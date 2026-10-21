#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
abggen.py -- 症例パラメータから「数式的に矛盾のない」動脈血液ガス分析を生成する
================================================================================
ecgen.py / labgen.py に続く第3のモダリティ。設計思想は共通：
  生成AIはパラメータJSONだけを書き、数値の算出と検証は決定論的コードが担う。

血液ガスが特別なのは、三つの主要項目が方程式で拘束されている点にある。
  pH = 6.1 + log10( HCO3⁻ / (0.03 × PaCO2) )     … Henderson-Hasselbalch
したがって pH を乱数で振ってはならない。本ツールでは
  「一次性障害 → 代償 → HCO3⁻ と PaCO2 が決まる → pH は計算で出る」
という順序を強制する。これにより pH 7.20 / HCO3 24 / PaCO2 40 のような
物理的にありえない血ガスが原理的に生成されない。

依存: numpy, matplotlib（図が不要なら numpy のみで動作）
使い方:
    python abggen.py case.json -o out    # out.html / out.png / out.json
    python abggen.py --demo              # 内蔵デモ症例を一括生成
    python abggen.py --profiles          # 病態プロファイル一覧
"""

import html as _html
import json
import math
import argparse
import numpy as np

# ============================================================================
# 基準範囲
# ============================================================================
ABG_REF = {
    "pH":    ("pH",          "",       7.35, 7.45, 3, "酸塩基"),
    "PaCO2": ("PaCO₂",       "mmHg",   35.0, 45.0, 1, "酸塩基"),
    "HCO3":  ("HCO₃⁻",       "mEq/L",  22.0, 26.0, 1, "酸塩基"),
    "BE":    ("Base Excess", "mEq/L",  -2.0,  2.0, 1, "酸塩基"),
    "PaO2":  ("PaO₂",        "mmHg",   80.0, 100.0, 1, "酸素化"),
    "SaO2":  ("SaO₂",        "%",      95.0, 99.0, 1, "酸素化"),
    "AaDO2": ("A-aDO₂",      "mmHg",    0.0, 15.0, 1, "酸素化"),
    "PF":    ("P/F比",       "",      400.0, 500.0, 0, "酸素化"),
    "Lac":   ("乳酸",        "mmol/L",  0.5,  2.0, 2, "代謝"),
    "AG":    ("アニオンギャップ", "mEq/L", 10.0, 14.0, 1, "代謝"),
    "cAG":   ("補正AG(Alb)", "mEq/L",  10.0, 14.0, 1, "代謝"),
    "Na":    ("Na⁺",         "mEq/L", 138.0, 145.0, 0, "電解質"),
    "K":     ("K⁺",          "mEq/L",   3.6,  4.8, 1, "電解質"),
    "Cl":    ("Cl⁻",         "mEq/L", 101.0, 108.0, 0, "電解質"),
    "Alb":   ("アルブミン",   "g/dL",    4.1,  5.1, 1, "電解質"),
    "Glu":   ("血糖",        "mg/dL",  73.0, 109.0, 0, "代謝"),
    "Hb":    ("ヘモグロビン", "g/dL",   13.7, 16.8, 1, "代謝"),
}

# ============================================================================
# 病態プロファイル
#   primary: 一次性障害
#   HCO3 / PaCO2 のうち「一次性で動く方」を指定し、他方は代償式で自動計算する。
#   compensation: "appropriate"(適切) / "acute" / "chronic" /
#                 "inadequate"(代償不全) / "none"(代償なし=混合性)
# ============================================================================
PROFILES = {
    "normal": {
        "_label": "正常", "primary": "none",
        "set": {}, "compensation": "appropriate",
    },
    "dka": {
        "_label": "糖尿病性ケトアシドーシス",
        "primary": "metabolic_acidosis", "compensation": "appropriate",
        "set": {"HCO3": 8.0, "Na": 134, "Cl": 94, "K": 5.6,
                "Glu": 486, "Lac": 1.8, "Alb": 4.2},
        "_ag_type": "high",
    },
    "lactic_acidosis": {
        "_label": "乳酸アシドーシス（敗血症性ショック）",
        "primary": "metabolic_acidosis", "compensation": "appropriate",
        "set": {"HCO3": 12.0, "Na": 138, "Cl": 104, "K": 4.4,
                "Lac": 8.6, "Alb": 2.6, "PaO2": 68.0, "FiO2": 0.4},
        "_ag_type": "high",
    },
    "diarrhea_acidosis": {
        "_label": "下痢によるAG正常性代謝性アシドーシス",
        "primary": "metabolic_acidosis", "compensation": "appropriate",
        "set": {"HCO3": 14.0, "Na": 138, "Cl": 116, "K": 2.9,
                "Lac": 1.2, "Alb": 3.8},
        "_ag_type": "normal",
    },
    "uremic_acidosis": {
        "_label": "尿毒症性アシドーシス（CKD G5）",
        "primary": "metabolic_acidosis", "compensation": "appropriate",
        "set": {"HCO3": 15.0, "Na": 136, "Cl": 104, "K": 6.1,
                "Lac": 1.4, "Alb": 3.2},
        "_ag_type": "high",
    },
    "vomiting_alkalosis": {
        "_label": "嘔吐による代謝性アルカローシス",
        "primary": "metabolic_alkalosis", "compensation": "appropriate",
        "set": {"HCO3": 38.0, "Na": 140, "Cl": 88, "K": 2.8,
                "Lac": 1.0, "Alb": 4.4},
        "_ag_type": "normal",
    },
    "copd_chronic": {
        "_label": "COPD 慢性呼吸性アシドーシス",
        "primary": "respiratory_acidosis", "compensation": "chronic",
        "set": {"PaCO2": 68.0, "Na": 140, "Cl": 96, "K": 4.2,
                "PaO2": 56.0, "FiO2": 0.21, "Alb": 3.9, "Hb": 17.8},
    },
    "copd_exacerbation": {
        "_label": "COPD増悪（慢性＋急性呼吸性アシドーシス）",
        "primary": "respiratory_acidosis", "compensation": "inadequate",
        "set": {"PaCO2": 88.0, "HCO3": 32.0, "Na": 139, "Cl": 94,
                "K": 4.6, "PaO2": 48.0, "FiO2": 0.28, "Alb": 3.6},
    },
    "opioid_overdose": {
        "_label": "急性呼吸性アシドーシス（薬剤性換気抑制）",
        "primary": "respiratory_acidosis", "compensation": "acute",
        "set": {"PaCO2": 74.0, "Na": 140, "Cl": 102, "K": 4.4,
                "PaO2": 54.0, "FiO2": 0.21, "Alb": 4.2},
    },
    "hyperventilation": {
        "_label": "急性呼吸性アルカローシス（過換気）",
        "primary": "respiratory_alkalosis", "compensation": "acute",
        "set": {"PaCO2": 24.0, "Na": 140, "Cl": 104, "K": 3.3,
                "PaO2": 108.0, "FiO2": 0.21, "Alb": 4.4},
    },
    "cirrhosis_alkalosis": {
        "_label": "肝硬変の慢性呼吸性アルカローシス",
        "primary": "respiratory_alkalosis", "compensation": "chronic",
        "set": {"PaCO2": 30.0, "Na": 129, "Cl": 94, "K": 3.4,
                "PaO2": 78.0, "FiO2": 0.21, "Alb": 2.4, "Lac": 2.4},
    },
    "pulmonary_embolism": {
        "_label": "肺塞栓症",
        "primary": "respiratory_alkalosis", "compensation": "acute",
        "set": {"PaCO2": 28.0, "Na": 139, "Cl": 103, "K": 3.8,
                "PaO2": 52.0, "FiO2": 0.21, "Alb": 4.0, "Lac": 2.2},
    },
    "ards": {
        "_label": "ARDS（重症呼吸不全）",
        "primary": "respiratory_acidosis", "compensation": "acute",
        "set": {"PaCO2": 52.0, "Na": 140, "Cl": 104, "K": 4.2,
                "PaO2": 62.0, "FiO2": 0.8, "Alb": 2.8, "Lac": 3.2},
    },
    "cardiac_arrest": {
        "_label": "心肺停止後（混合性アシドーシス）",
        "primary": "metabolic_acidosis", "compensation": "none",
        "set": {"HCO3": 9.0, "PaCO2": 62.0, "Na": 138, "Cl": 100,
                "K": 6.4, "Lac": 14.2, "PaO2": 58.0, "FiO2": 1.0, "Alb": 3.0},
        "_ag_type": "high",
    },
    "severe_pancreatitis": {
        "_label": "重症急性膵炎（乳酸アシドーシス＋呼吸性アルカローシス）",
        "primary": "metabolic_acidosis", "compensation": "none",
        "set": {"HCO3": 16.0, "PaCO2": 26.0, "Na": 136, "Cl": 100,
                "K": 3.6, "Lac": 5.4, "PaO2": 66.0, "FiO2": 0.35, "Alb": 2.4},
        "_ag_type": "high",
    },
    "vomiting_plus_aki": {
        "_label": "幽門狭窄（代謝性アルカローシス＋AG開大性アシドーシス）",
        "primary": "metabolic_alkalosis", "compensation": "appropriate",
        "set": {"HCO3": 30.0, "Na": 137, "Cl": 82, "K": 2.6,
                "Lac": 3.8, "Alb": 3.4},
        "_ag_type": "high",
    },
}


def default_params():
    return {
        "patient": {"id": "ABG-001", "age": 60, "sex": "M", "note": ""},
        "profile": "normal",
        "overrides": {},              # {"PaCO2": 55, "Lac": 4.2, "FiO2": 0.5}
        "sample": "動脈血",
        "temperature_c": 37.0,
        "altitude_patm": 760.0,
        "seed": 42,
        "variability": 0.04,
        "diagnosis_label": "",
        "show_interpretation": True,  # 演習で伏せるなら false
        "collected_at": "",
    }


# ============================================================================
# 生理学の式
# ============================================================================
def ph_from(hco3, paco2):
    """Henderson-Hasselbalch。pHは常にここでのみ決まる。"""
    return 6.1 + math.log10(max(hco3, 0.1) / (0.03 * max(paco2, 1.0)))


def base_excess(hco3, ph):
    """Van Slyke式による標準塩基過剰"""
    return 0.9287 * (hco3 - 24.4 + 14.83 * (ph - 7.4))


def sao2_from_pao2(pao2):
    """Severinghaus式による酸素飽和度の推定"""
    p = max(pao2, 1.0)
    return 100.0 / (23400.0 / (p ** 3 + 150.0 * p) + 1.0)


def aado2(pao2, paco2, fio2, patm=760.0, rq=0.8):
    """肺胞気-動脈血酸素分圧較差"""
    pao2_alv = fio2 * (patm - 47.0) - paco2 / rq
    return pao2_alv - pao2


def expected_aado2(age):
    """年齢による予測A-aDO2（室内気）"""
    return age / 4.0 + 4.0


def expected_compensation(primary, hco3, paco2, mode):
    """
    一次性障害に対する予測代償値と許容幅を返す。
    戻り値: (対象項目, 予測値, 許容幅, 用いた式の名称)
    """
    if primary == "metabolic_acidosis":
        # Winterの式
        return ("PaCO2", 1.5 * hco3 + 8.0, 2.0, "Winterの式: PaCO₂ = 1.5×HCO₃⁻ + 8 ± 2")
    if primary == "metabolic_alkalosis":
        return ("PaCO2", 0.7 * hco3 + 20.0, 5.0, "PaCO₂ = 0.7×HCO₃⁻ + 20 ± 5")
    if primary == "respiratory_acidosis":
        d = (paco2 - 40.0) / 10.0
        if mode == "chronic":
            return ("HCO3", 24.0 + 3.5 * d, 2.0,
                    "慢性呼吸性アシドーシス: ΔHCO₃⁻ = 3.5 × ΔPaCO₂/10")
        return ("HCO3", 24.0 + 1.0 * d, 2.0,
                "急性呼吸性アシドーシス: ΔHCO₃⁻ = 1.0 × ΔPaCO₂/10")
    if primary == "respiratory_alkalosis":
        d = (40.0 - paco2) / 10.0
        if mode == "chronic":
            return ("HCO3", 24.0 - 4.0 * d, 2.0,
                    "慢性呼吸性アルカローシス: ΔHCO₃⁻ = 4.0 × ΔPaCO₂/10")
        return ("HCO3", 24.0 - 2.0 * d, 2.0,
                "急性呼吸性アルカローシス: ΔHCO₃⁻ = 2.0 × ΔPaCO₂/10")
    return (None, None, None, None)


# ============================================================================
# 生成本体
# ============================================================================
def generate_values(p):
    rng = np.random.default_rng(p.get("seed", 42))
    prof = PROFILES.get(p.get("profile", "normal"), PROFILES["normal"])
    var = float(p.get("variability", 0.04))
    age = int(p["patient"].get("age", 60))
    patm = float(p.get("altitude_patm", 760.0))

    # --- 既定値（正常）---
    v = {"HCO3": 24.0, "PaCO2": 40.0, "Na": 140.0, "K": 4.2, "Cl": 104.0,
         "Alb": 4.4, "Lac": 1.0, "Glu": 92.0, "Hb": 14.6,
         "PaO2": 92.0, "FiO2": 0.21}
    for k in ("Na", "Cl", "K", "Glu", "Hb"):
        v[k] += rng.normal(0, abs(v[k]) * var * 0.3)

    # --- プロファイル適用 ---
    st = dict(prof.get("set", {}))
    for k, val in st.items():
        v[k] = float(val)

    primary = prof.get("primary", "none")
    comp = prof.get("compensation", "appropriate")

    # --- 代償の計算：一次性で動かない方の値を式から決める ---
    comp_note = None
    if primary != "none" and comp != "none":
        target, pred, tol, formula = expected_compensation(
            primary, v["HCO3"], v["PaCO2"], comp)
        if target and target not in st:
            if comp == "inadequate":
                # 代償不全：予測値の6割程度しか動かない
                base = 40.0 if target == "PaCO2" else 24.0
                v[target] = base + (pred - base) * 0.6
            else:
                v[target] = pred + rng.normal(0, tol * 0.4)
        comp_note = formula

    # --- PaO2 の物理的結合 ---
    #   プロファイルがPaO2を明示していない場合、既定値を使うとPaCO2の変化と
    #   矛盾する（過換気なのにPaO2が正常のままA-aDO2だけ開大する等）。
    #   肺胞気式で PAO2 を求め、年齢相当のA-aDO2を引いて PaO2 を導出する。
    if "PaO2" not in st and "PaO2" not in (p.get("overrides") or {}):
        pao2_alv = v["FiO2"] * (patm - 47.0) - v["PaCO2"] / 0.8
        v["PaO2"] = pao2_alv - expected_aado2(age)

    # --- ユーザー上書き ---
    for k, val in (p.get("overrides") or {}).items():
        v[k] = float(val)

    # --- 派生値（★乱数を使わない）---
    v["pH"] = ph_from(v["HCO3"], v["PaCO2"])
    v["BE"] = base_excess(v["HCO3"], v["pH"])
    v["SaO2"] = sao2_from_pao2(v["PaO2"])
    v["AaDO2"] = aado2(v["PaO2"], v["PaCO2"], v["FiO2"], patm)
    v["PF"] = v["PaO2"] / max(v["FiO2"], 0.01)
    v["AG"] = v["Na"] - v["Cl"] - v["HCO3"]
    v["cAG"] = v["AG"] + 2.5 * (4.0 - v["Alb"])
    return v, prof, comp_note


def analyze(v, p, prof, comp_note):
    """判読ステップを自動生成する（教育用の解答例に相当）"""
    age = int(p["patient"].get("age", 60))
    steps = []
    ph, hco3, paco2 = v["pH"], v["HCO3"], v["PaCO2"]

    # Step 1 アシデミア/アルカレミア
    if ph < 7.35:
        s1 = f"pH {ph:.3f} → アシデミア"
    elif ph > 7.45:
        s1 = f"pH {ph:.3f} → アルカレミア"
    else:
        s1 = f"pH {ph:.3f} → 正常範囲（代償済み、あるいは混合性障害の可能性）"
    steps.append(("① pHの評価", s1))

    # Step 2 一次性障害
    if hco3 < 22 and paco2 < 35:
        s2 = "HCO₃⁻低下・PaCO₂低下 → 代謝性アシドーシス（呼吸性代償）が主体"
    elif hco3 > 26 and paco2 > 45:
        s2 = "HCO₃⁻上昇・PaCO₂上昇 → 代謝性アルカローシスまたは呼吸性アシドーシス"
    elif paco2 > 45:
        s2 = "PaCO₂上昇 → 呼吸性アシドーシス"
    elif paco2 < 35:
        s2 = "PaCO₂低下 → 呼吸性アルカローシス"
    elif hco3 < 22:
        s2 = "HCO₃⁻低下 → 代謝性アシドーシス"
    elif hco3 > 26:
        s2 = "HCO₃⁻上昇 → 代謝性アルカローシス"
    else:
        s2 = "HCO₃⁻・PaCO₂ともに正常範囲"
    steps.append(("② 一次性障害の同定",
                  s2 + f"（HCO₃⁻ {hco3:.1f}, PaCO₂ {paco2:.1f}）"))

    # Step 3 代償の妥当性
    primary = prof.get("primary", "none")
    if primary != "none":
        target, pred, tol, formula = expected_compensation(
            primary, hco3, paco2, prof.get("compensation", "appropriate"))
        if target:
            actual = v[target]
            unit = "mmHg" if target == "PaCO2" else "mEq/L"
            lo, hi = pred - tol, pred + tol
            if lo <= actual <= hi:
                verdict = "予測範囲内 → 単一の障害として矛盾しない"
            elif actual < lo:
                verdict = ("予測より低値 → " +
                           ("呼吸性アルカローシスの合併" if target == "PaCO2"
                            else "代謝性アシドーシスの合併"))
            else:
                verdict = ("予測より高値 → " +
                           ("呼吸性アシドーシスの合併" if target == "PaCO2"
                            else "代謝性アルカローシスの合併"))
            steps.append(("③ 代償の評価",
                          f"{formula}<br>予測 {pred:.1f} ({lo:.1f}–{hi:.1f}) "
                          f"{unit} / 実測 {actual:.1f} {unit} → {verdict}"))

    # Step 4 アニオンギャップ
    ag, cag = v["AG"], v["cAG"]
    s4 = f"AG {ag:.1f} mEq/L"
    if v["Alb"] < 4.0:
        s4 += f"／Alb {v['Alb']:.1f} g/dL による補正AG {cag:.1f} mEq/L"
        s4 += "（低Alb血症ではAGが過小評価されるため補正が必須）"
    if cag > 16:
        s4 += " → AG開大性代謝性アシドーシスあり"
    elif cag < 10:
        s4 += " → AG低下（低Alb・多発性骨髄腫など）"
    else:
        s4 += " → AGは正常範囲"
    steps.append(("④ アニオンギャップ", s4))

    # Step 5 デルタ比（AG開大時のみ）
    if cag > 16 and hco3 < 24:
        delta_ratio = (cag - 12.0) / (24.0 - hco3)
        if delta_ratio < 0.4:
            d = "AG正常性代謝性アシドーシスの合併"
        elif delta_ratio < 1.0:
            d = "AG開大性＋AG正常性アシドーシスの混在"
        elif delta_ratio <= 2.0:
            d = "純粋なAG開大性代謝性アシドーシスとして矛盾しない"
        else:
            d = "代謝性アルカローシスまたは慢性呼吸性アシドーシスの合併"
        steps.append(("⑤ デルタ比（ΔAG/ΔHCO₃⁻）",
                      f"補正AGを用いて {delta_ratio:.2f} → {d}"))
    elif cag > 16 and hco3 >= 24:
        steps.append(("⑤ デルタ比の考え方",
                      f"補正AGが {cag:.1f} と開大しているにもかかわらず "
                      f"HCO₃⁻ {hco3:.1f} mEq/L が低下していない。"
                      "AG開大分だけHCO₃⁻は下がるはずであり、"
                      "代謝性アルカローシスの合併を強く示唆する"))

    # Step 6 酸素化
    exp_aa = expected_aado2(age)
    s6 = (f"PaO₂ {v['PaO2']:.0f} mmHg（FiO₂ {v['FiO2']:.2f}）／"
          f"A-aDO₂ {v['AaDO2']:.1f} mmHg（年齢予測 {exp_aa:.1f}）")
    if v["AaDO2"] > exp_aa + 10:
        s6 += " → A-aDO₂開大：換気血流不均衡・シャント・拡散障害を示唆"
    elif v["PaCO2"] > 45:
        s6 += " → A-aDO₂は開大なし：低酸素血症は肺胞低換気で説明できる"
    elif v["PaO2"] < 80:
        s6 += " → A-aDO₂は開大なし（吸入酸素分圧の低下等を検討）"
    else:
        s6 += " → A-aDO₂・PaO₂ともに問題なし"
    if v["FiO2"] > 0.21:
        s6 += f"／P/F比 {v['PF']:.0f}"
        if v["PF"] < 100:
            s6 += "（重症ARDS域）"
        elif v["PF"] < 200:
            s6 += "（中等症ARDS域）"
        elif v["PF"] < 300:
            s6 += "（軽症ARDS域）"
    steps.append(("⑥ 酸素化の評価", s6))

    # Step 7 乳酸
    if v["Lac"] > 2.0:
        sev = "著明高値（組織低灌流を強く示唆）" if v["Lac"] >= 4.0 else "軽度高値"
        steps.append(("⑦ 乳酸", f"{v['Lac']:.2f} mmol/L → {sev}"))
    return steps


# ============================================================================
# 整合性検証
# ============================================================================
def validate(v, p, prof):
    w = []
    # Henderson-Hasselbalch の再検算
    ph_calc = ph_from(v["HCO3"], v["PaCO2"])
    if abs(ph_calc - v["pH"]) > 0.005:
        w.append(f"pHがHenderson-Hasselbalch式と不一致"
                 f"（実測 {v['pH']:.3f} vs 計算 {ph_calc:.3f}）")
    if not (6.60 <= v["pH"] <= 7.90):
        w.append(f"pH {v['pH']:.3f} は生存とほぼ両立しない値")
    if v["HCO3"] <= 0 or v["PaCO2"] <= 0:
        w.append("HCO₃⁻またはPaCO₂が非正値")
    if v["PaO2"] > v["FiO2"] * (760 - 47):
        w.append("PaO₂が吸入酸素分圧の理論上限を超えている")
    if v["AaDO2"] < -5:
        w.append("A-aDO₂が負値（PaO₂とFiO₂の組み合わせが不整合）")
    if v["AG"] < 3:
        w.append(f"アニオンギャップ {v['AG']:.1f} は低すぎる（Na/Cl/HCO₃⁻を確認）")
    if v["Lac"] > 4.0 and v["cAG"] < 14:
        w.append("著明な乳酸高値に対しAGが開大していない（Na/Clを確認）")
    if v["K"] > 6.0 and v["pH"] > 7.45:
        w.append("アルカレミアと高K血症の併存は要根拠")
    if v["pH"] < 7.20 and v["HCO3"] > 20 and v["PaCO2"] < 50:
        w.append("重度アシデミアに対しHCO₃⁻・PaCO₂の変化が乏しい")
    sao2_calc = sao2_from_pao2(v["PaO2"])
    if abs(sao2_calc - v["SaO2"]) > 1.0:
        w.append("SaO₂がPaO₂と不整合")
    return w


# ============================================================================
# 出力
# ============================================================================
def flag(val, lo, hi):
    if val > hi:
        return "H"
    if val < lo:
        return "L"
    return ""


def render_html(v, p, steps, warns, prof):
    pt = p["patient"]
    sex = "男性" if pt.get("sex") == "M" else "女性"
    css = """
    body{font-family:"Hiragino Sans","Noto Sans JP","Yu Gothic",sans-serif;
         margin:24px;color:#111;background:#fff;max-width:820px}
    h1{font-size:17px;border-bottom:2px solid #333;padding-bottom:6px;margin:0 0 4px}
    .meta{font-size:12px;color:#444;margin-bottom:14px}
    h2{font-size:13px;background:#eef2f6;padding:4px 8px;margin:16px 0 2px;
       border-left:4px solid #4a6785}
    table{border-collapse:collapse;width:100%;font-size:12px}
    td,th{padding:3px 8px;border-bottom:1px solid #e3e3e3;text-align:right}
    th{background:#fafafa;font-size:11px;color:#555}
    td.n,th.n{text-align:left;width:34%}
    td.v{font-weight:600;width:14%}td.f{width:6%;font-weight:700}
    .H{color:#c0392b}.L{color:#1a5f9e}
    .unit{color:#666;font-size:11px;width:16%;text-align:left}
    .ref{color:#666;font-size:11px;width:24%}
    .step{font-size:12.5px;margin:6px 0;padding:7px 12px;background:#f7f9fb;
          border-left:3px solid #7d9bb5;line-height:1.6}
    .step b{color:#2c4a63}
    .warn{font-size:12px;background:#fff6e5;border-left:4px solid #e0a800;
          padding:8px 12px;margin-top:12px}
    .ok{font-size:12px;color:#2d7a3e;margin-top:12px}
    .note{font-size:11px;color:#777;margin-top:18px;border-top:1px solid #ddd;
          padding-top:8px;line-height:1.6}
    .hh{font-size:12px;background:#f2f7f2;border-left:4px solid #6a9b6a;
        padding:8px 12px;margin-top:8px;font-family:monospace}
    """
    # 症例JSON（生成AIが作ることもある）の文字列は、HTML として解釈されないようエスケープする
    e = lambda x: _html.escape(str(x))
    h = ["<!DOCTYPE html><html lang='ja'><head><meta charset='utf-8'>",
         # スクリプトや外部読み込みを一切許可しない（報告書は表示だけ）
         "<meta http-equiv='Content-Security-Policy' content=\"default-src 'none'; style-src 'unsafe-inline'; img-src data:\">",
         f"<title>血液ガス分析 {e(pt.get('id',''))}</title><style>{css}</style></head><body>",
         "<h1>動脈血液ガス分析報告書</h1>"]
    meta = f"ID: {e(pt.get('id',''))}　{e(pt.get('age',''))}歳　{e(sex)}　検体: {e(p.get('sample','動脈血'))}"
    if p.get("collected_at"):
        meta += f"　採取: {e(p['collected_at'])}"
    if pt.get("note"):
        meta += f"　／　{e(pt['note'])}"
    h.append(f"<div class='meta'>{meta}</div>")

    for panel in ["酸塩基", "酸素化", "代謝", "電解質"]:
        rows = [(k, r) for k, r in ABG_REF.items() if r[5] == panel]
        if not rows:
            continue
        h.append(f"<h2>{panel}</h2><table>")
        h.append("<tr><th class='n'>項目</th><th>結果</th><th></th>"
                 "<th class='unit'>単位</th><th class='ref'>基準範囲</th></tr>")
        for k, (name, unit, lo, hi, dec, _) in rows:
            if k not in v:
                continue
            val = v[k]
            fl = flag(val, lo, hi)
            h.append(f"<tr><td class='n'>{name}</td>"
                     f"<td class='v {fl}'>{val:.{dec}f}</td>"
                     f"<td class='f {fl}'>{fl}</td><td class='unit'>{unit}</td>"
                     f"<td class='ref'>{lo:.{dec}f} – {hi:.{dec}f}</td></tr>")
        h.append("</table>")

    h.append(f"<div class='hh'>Henderson-Hasselbalch 検算：<br>"
             f"pH = 6.1 + log₁₀( {v['HCO3']:.1f} / (0.03 × {v['PaCO2']:.1f}) ) "
             f"= {ph_from(v['HCO3'], v['PaCO2']):.3f}　✓</div>")

    if p.get("show_interpretation", True):
        h.append("<h2>判読ステップ</h2>")
        for title, body in steps:
            h.append(f"<div class='step'><b>{title}</b><br>{body}</div>")
        if p.get("diagnosis_label") or prof.get("_label"):
            lab = p.get("diagnosis_label") or prof.get("_label")
            h.append(f"<div class='note'><b>想定病態：</b>{e(lab)}</div>")
    else:
        h.append("<div class='note'>※ 判読は演習用に非表示にしています。</div>")

    if warns:
        h.append("<div class='warn'><b>整合性チェック：要確認</b><ul>"
                 + "".join(f"<li>{x}</li>" for x in warns) + "</ul></div>")
    else:
        h.append("<div class='ok'>✓ 整合性チェック：数式上の矛盾は検出されませんでした</div>")

    h.append("<div class='note'>本データは教育目的で生成された合成検査値であり、"
             "実在の患者に由来するものではありません。代償予測式・基準範囲は"
             "成書に基づく一般値であり、運用時は自施設の方針に合わせてください。</div>")
    h.append("</body></html>")
    return "\n".join(h)


def plot_davenport(v, p, outpath):
    """Davenport図（pH–HCO₃⁻ 平面）に症例をプロットする教育用の図"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib import font_manager
        for f in ["Noto Sans CJK JP", "Noto Sans JP", "IPAexGothic",
                  "Hiragino Sans", "Yu Gothic", "MS Gothic"]:
            if any(f == x.name for x in font_manager.fontManager.ttflist):
                plt.rcParams["font.family"] = f
                break
    except Exception:
        return None

    fig, ax = plt.subplots(figsize=(6.4, 4.8), dpi=160)
    ph = np.linspace(6.9, 7.75, 300)
    for pco2 in [10, 20, 30, 40, 50, 60, 80, 100]:
        hco3 = 0.03 * pco2 * (10 ** (ph - 6.1))
        ax.plot(ph, hco3, lw=0.8, color="#9bb7cc")
        if hco3[-1] >= 44:                      # 上端で曲線が抜ける場合
            idx = int(np.argmin(np.abs(hco3 - 44)))
            ax.text(ph[idx], 45.0, f"{pco2}", fontsize=7, color="#6f8ea6",
                    ha="center")
        else:                                    # 右端で抜ける場合
            ax.text(7.755, hco3[-1], f"{pco2}", fontsize=7, color="#6f8ea6",
                    va="center", ha="left")
    ax.axhspan(22, 26, color="#dfe9c8", alpha=0.5, zorder=0)
    ax.axvspan(7.35, 7.45, color="#dfe9c8", alpha=0.5, zorder=0)
    ax.plot(v["pH"], v["HCO3"], "o", ms=11, color="#c0392b",
            markeredgecolor="white", markeredgewidth=1.4, zorder=5)
    ax.annotate(f"pH {v['pH']:.2f}\nHCO3- {v['HCO3']:.1f}\nPaCO2 {v['PaCO2']:.0f}",
                (v["pH"], v["HCO3"]), textcoords="offset points",
                xytext=(12, 10), fontsize=9,
                bbox=dict(boxstyle="round,pad=0.35", fc="#fff4f2",
                          ec="#c0392b", lw=0.8))
    ax.set_xlim(6.9, 7.75)
    ax.set_ylim(0, 48)
    ax.set_xlabel("pH")
    ax.set_ylabel("HCO3- (mEq/L)")
    ax.set_title("Davenport図（等PaCO2線は mmHg）", fontsize=11)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)
    return outpath


def generate(case, outbase):
    p = default_params()
    for k, val in case.items():
        if isinstance(val, dict) and isinstance(p.get(k), dict):
            p[k].update(val)
        else:
            p[k] = val

    v, prof, comp_note = generate_values(p)
    if not p.get("diagnosis_label"):
        p["diagnosis_label"] = prof.get("_label", "")
    steps = analyze(v, p, prof, comp_note)
    warns = validate(v, p, prof)

    with open(f"{outbase}.html", "w", encoding="utf-8") as f:
        f.write(render_html(v, p, steps, warns, prof))
    png = plot_davenport(v, p, f"{outbase}.png")

    out = {"params": p,
           "values": {k: round(val, 4) for k, val in v.items()},
           "interpretation": [{"step": t, "body": b} for t, b in steps],
           "validation_warnings": warns}
    with open(f"{outbase}.json", "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=2)
    return f"{outbase}.html", png, warns


DEMO = {
    "G01_normal":        {"profile": "normal",
                          "patient": {"id": "ABG-01", "age": 40, "sex": "M"}},
    "G02_dka":           {"profile": "dka",
                          "patient": {"id": "ABG-02", "age": 24, "sex": "F",
                                      "note": "意識障害・Kussmaul呼吸"}},
    "G03_septic_lactic": {"profile": "lactic_acidosis",
                          "patient": {"id": "ABG-03", "age": 72, "sex": "M",
                                      "note": "敗血症性ショック"}},
    "G04_diarrhea":      {"profile": "diarrhea_acidosis",
                          "patient": {"id": "ABG-04", "age": 33, "sex": "F",
                                      "note": "1週間続く水様下痢"}},
    "G05_vomiting":      {"profile": "vomiting_alkalosis",
                          "patient": {"id": "ABG-05", "age": 58, "sex": "M",
                                      "note": "反復する嘔吐"}},
    "G06_copd":          {"profile": "copd_chronic",
                          "patient": {"id": "ABG-06", "age": 76, "sex": "M",
                                      "note": "安定期COPD"}},
    "G07_copd_exac":     {"profile": "copd_exacerbation",
                          "patient": {"id": "ABG-07", "age": 78, "sex": "M",
                                      "note": "CO₂ナルコーシス"}},
    "G08_cirrhosis":     {"profile": "cirrhosis_alkalosis",
                          "patient": {"id": "LAB-03", "age": 62, "sex": "M",
                                      "note": "labgen.py L03_cirrhosis と対応"}},
    "G09_pancreatitis":  {"profile": "severe_pancreatitis",
                          "patient": {"id": "LAB-05", "age": 51, "sex": "M",
                                      "note": "labgen.py L05 と対応"}},
    "G10_arrest":        {"profile": "cardiac_arrest",
                          "patient": {"id": "ABG-10", "age": 68, "sex": "M",
                                      "note": "ROSC直後"}},
    "G11_pe":            {"profile": "pulmonary_embolism",
                          "patient": {"id": "ABG-11", "age": 64, "sex": "F",
                                      "note": "突然の呼吸困難"}},
    "G12_pyloric":       {"profile": "vomiting_plus_aki",
                          "patient": {"id": "ABG-12", "age": 70, "sex": "M",
                                      "note": "幽門狭窄・脱水"}},
}


def main():
    ap = argparse.ArgumentParser(description="症例パラメータから動脈血液ガスを生成")
    ap.add_argument("case", nargs="?")
    ap.add_argument("-o", "--out", default="abg")
    ap.add_argument("--demo", action="store_true")
    ap.add_argument("--schema", action="store_true")
    ap.add_argument("--profiles", action="store_true")
    a = ap.parse_args()

    if a.profiles:
        for k, val in PROFILES.items():
            print(f"{k:24s} {val.get('_label','')}")
        return
    if a.schema:
        print(json.dumps(default_params(), ensure_ascii=False, indent=2))
        return
    if a.demo:
        for name, case in DEMO.items():
            html, png, w = generate(case, name)
            print(("⚠" if w else "✓") + f" {html}" +
                  (f"  ({len(w)}件の要確認)" if w else ""))
        return
    if not a.case:
        ap.error("症例JSONを指定するか --demo / --schema / --profiles を使ってください")
    with open(a.case, encoding="utf-8") as f:
        case = json.load(f)
    html, png, w = generate(case, a.out)
    print("generated:", html, png)
    for x in w:
        print("  ⚠", x)


if __name__ == "__main__":
    main()
