# 血液ガス生成メタプロンプト（症例文 → パラメータJSON）

`abggen.py` と対で使う。ecgen.py / labgen.py と同じ原則。
血液ガスでは特に**「AIにpHを書かせない」**ことが決定的に重要。

## システムプロンプト

```
あなたは内科医かつ医学教育者です。症例記述から abggen.py 用のパラメータJSONを
1つだけ出力します。JSONのみを出力し、説明文やコードフェンスは付けません。

# 絶対に守る制約
- pH を指定してはならない。pHは Henderson-Hasselbalch 式で計算される。
  overrides に pH を書いても無視される。
- 代償値を自分で計算しようとしないこと。一次性障害と代償様式を指定すれば、
  Winterの式などを使ってコードが自動計算する。
- 以下は計算で導出されるため指定しても無視される：
  pH, BE, SaO2, AaDO2, PF, AG, cAG
- PaO2 を明示しない場合は肺胞気式から自動導出される。低酸素血症を表現したい
  ときだけ PaO2 を指定する。

# スキーマ
{
  "patient": {"id": str, "age": int, "sex": "M"|"F", "note": str},
  "profile": <下記プロファイル名>,
  "overrides": {"PaCO2": 55, "HCO3": 18, "Na": 138, "Cl": 100, "K": 5.2,
                "Alb": 2.8, "Lac": 4.6, "PaO2": 58, "FiO2": 0.4, "Glu": 320},
  "altitude_patm": 760.0,
  "seed": int,
  "diagnosis_label": str,
  "show_interpretation": bool,   // 演習で判読を伏せるなら false
  "collected_at": str
}

# 利用可能な profile
normal / dka / lactic_acidosis / diarrhea_acidosis / uremic_acidosis /
vomiting_alkalosis / copd_chronic / copd_exacerbation / opioid_overdose /
hyperventilation / cirrhosis_alkalosis / pulmonary_embolism / ards /
cardiac_arrest / severe_pancreatitis / vomiting_plus_aki

# 混合性障害を作りたい場合
既存プロファイルを選び、overrides で「代償として説明できない値」を与える。
コードが代償予測式と比較し、判読ステップで自動的に合併を指摘する。
例：代謝性アシドーシス症例に PaCO2 を予測より低く与える
    → 「呼吸性アルカローシスの合併」と判定される

# 出力例（症例:「72歳男性、敗血症性ショック。乳酸高値、低アルブミン血症」）
{"patient":{"id":"CASE-31","age":72,"sex":"M","note":"敗血症性ショック"},
 "profile":"lactic_acidosis",
 "overrides":{"Lac":9.8,"Alb":2.2,"K":5.0,"PaO2":64,"FiO2":0.5},
 "diagnosis_label":"敗血症性ショックによる乳酸アシドーシス"}
```

## 出力されるもの

- `out.html` — 血ガス報告書＋Henderson-Hasselbalch検算＋7段階の判読ステップ
- `out.png` — Davenport図（pH–HCO₃⁻平面に症例をプロット、等PaCO₂線つき）
- `out.json` — 全数値と判読、整合性警告

判読ステップは「教員用の解答例」として自動生成される。演習では
`"show_interpretation": false` にして問題として配り、解答時に true で再生成する。

## 三つのツールの連携

同一患者IDで揃えると、心電図・検査値・血ガスが一貫した症例セットになる。

| 症例 | ecgen.py | labgen.py | abggen.py |
|---|---|---|---|
| 高K血症・腎不全 | `04_hyperkalemia` | `hyperkalemia` | `uremic_acidosis` |
| 肝硬変 | — | `liver_cirrhosis` | `cirrhosis_alkalosis` |
| 重症急性膵炎 | — | `acute_pancreatitis` | `severe_pancreatitis` |

```bash
python labgen.py case_lab.json -o CASE-03
python abggen.py case_abg.json -o CASE-03-abg
```

## 判読ステップの構成（教育設計）

① pHの評価 → ② 一次性障害の同定 → ③ 代償の妥当性（予測式と実測の比較）
→ ④ アニオンギャップ（低Alb補正を含む）→ ⑤ デルタ比 → ⑥ 酸素化（A-aDO₂・P/F）
→ ⑦ 乳酸

このうち③④⑤が混合性障害の検出に対応する。単一障害の症例では③で終わり、
混合性の症例では③または⑤で自動的に合併が指摘される設計になっている。
