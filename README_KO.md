# 2D 급수축 정상 유동

입자가 쌓이기 전의 채널에서 Gly75, PEO Cross, 2모드 지수형 PTT 유동을 구합니다. 단차의 바깥 코너에서 목 끝으로 향하는 거리 `s`에 대해 2 µm 측정선의 벽 방향 속도와 같은 위치의 벽 접선응력을 비교합니다.

## 설치와 첫 계산

GitHub `main`의 기준 커밋은 `df22a9f6bd6607754c1fa10ade7272133e6f71b6`입니다. 압축 파일의 경로는 저장소 루트 기준입니다. Gautschi에서 기존 프로젝트 위치와 Foundation OpenFOAM-8 설정을 사용합니다. 기존 작업이 같은 실행 파일을 사용 중이면 해당 작업이 끝난 뒤 빌드합니다.

```bash
cd /home/lee6456/Deposition
unzip -o Test_Deposit_steady2d_overlay.zip
build_id=$(sbatch --parsable build_dep.sbatch)
bash submit_comparison.sh 0.5 coarse --after-build "${build_id%%;*}"
```

기본 Slurm 계정은 현재 Git의 `lilly-nan6`입니다. 로컬에서 다른 allocation을 사용하면 압축을 푼 뒤 `build_dep.sbatch`, `run_dep.sbatch`, `postprocess_dep.sbatch`의 계정을 맞춥니다. 빌드에는 기존의 고정 버전 rheoTool 소스와 Eigen을 사용합니다. 그림을 만드는 Python 환경에는 `requirements-plot.txt`가 필요합니다.

0.5 mL/h에서는 세 모델을 제출하고 모두 수렴하면 후처리를 시작합니다. 결과 그림은 `comparison/steady2d/Q0p5_coarse/step_wall_comparison.png` 및 같은 이름의 PDF입니다. 다음 두 명령은 각각 사진과 같은 0.3 mL/h 계산과 0.5 mL/h 격자 비교입니다. 동일 설정의 수렴한 기본 격자 결과를 재사용합니다.

```bash
bash submit_comparison.sh 0.3 coarse
bash submit_comparison.sh 0.5 coarse --compare-mesh medium
```

격자 비교의 결과 폴더는 `comparison/steady2d/Q0p5_coarse_vs_medium/`입니다. 이후 `bash submit_comparison.sh 0.7 coarse`로 같은 절차를 실행할 수 있습니다. 단일 모델은 `sbatch run_dep.sbatch peo 0.5 coarse`로 실행합니다. 중단되거나 반복 상한에 도달한 계산은 `sbatch run_dep.sbatch peo 0.5 coarse --resume --max-iters 5000`으로 마지막 완전한 체크포인트부터 추가 반복합니다. 각 단계의 로그와 `postProcessing/deposition/metrics.csv`에서 미수렴 원인을 확인할 수 있습니다.

## 형상과 방정식

- 흐름 방향은 +x이고 단차는 x=0입니다. 상·하류 길이는 각각 900 µm, 폭은 각각 180/45 µm입니다. 위쪽 반 영역만 사용하여 중심선 y=0에 대칭 조건을 둡니다. 각 벽의 속도는 0입니다. 깊이 방향 한 셀의 두 면은 `empty`입니다.
- `height=100 µm`는 실험 유량 환산에 사용합니다. 단위 깊이당 반쪽 채널 유량은 `q=Q/(2h)`이고 0.5 mL/h에서 6.9444444×10⁻⁷ m²/s입니다. 격자의 수치적 두께는 1 µm입니다. 계산된 면 유량을 이 두께로 나누어 선 유량을 얻습니다.
- 평면 포물선 입구 속도를 각 입구 면 구간에서 평균하여 이산 유량을 맞춥니다. 각 ePTT 모드의 입구 응력과 로그 conformation은 주어진 포물선의 국소 전단률에 대응하는 정상 단순 전단 해에서 계산합니다. 출구에는 전체 Cauchy traction이 0인 조건을 적용하며 운동학적 압력과 속도의 법선구배를 함께 갱신합니다.
- 운동량 식은 정상 Stokes 식입니다. ePTT는 정상 상태에서도 로그 conformation의 공간 대류, 유동에 의한 변형, 이완을 포함합니다. Picard 반복에서 대각 안정화 항을 양변에 더하며, 수렴한 방정식에서 그 기여는 0입니다. 출력 디렉터리 번호는 반복 번호입니다.
- 기본 물성은 `config/experiment.json`에 있습니다. Gly75의 점도는 0.035 Pa·s입니다. PEO의 두 이완시간은 40.0000/23.4887 ms이고 용매 점도는 두 모드 합계로 0.00695 Pa·s입니다. Cross는 논문 식을 사용합니다. Cross와 ePTT의 **전단점도 적합 오차**는 유량·위치별 탄성 효과 해석에 함께 고려해야 합니다. 예를 들어 1372 s⁻¹에서 각각 0.03102, 0.03073 Pa·s이고 10 s⁻¹에서는 0.16716, 0.18853 Pa·s입니다.

## 벽과 와류 결과 읽기

| 출력 | 정의와 용도 |
| --- | --- |
| `step_wall_comparison.png/.pdf` | 위 패널은 단차 앞 x=−2 µm에서 +x 방향 속도, 아래 패널은 같은 s의 부호 있는 접선응력입니다. s=0은 바깥 코너, 67.5 µm는 목 끝입니다. |
| `wall_profiles.csv` | 세 벽에 대한 세 모델의 벽 방향 속도·벽 접선응력과 PEO 용매·모드별 응력입니다. 표면에서 밖으로 향하는 유체영역 법선 n을 쓰며 접선 t는 상·하류 +x, 단차 −y입니다. 출력은 t·(τ·n)입니다. 전단 크기 비교에는 절댓값을 사용합니다. |
| `approach_intervals.csv` | 벽 쪽 속도가 양수인 연속 구간의 양 끝, 폭, 속도 최대 위치, 그 위치의 응력, 구간별 ∑max(u·n,0)Δs 및 q에 대한 비율입니다. 구간 끝은 해당 격자 면의 경계로 정합니다. |
| `wall_integrals.csv`, `straight_wall_check.csv` | 전체 벽과 수축부 인접 직선벽(단차에서 180 µm 이내)의 값을 구분합니다. 입·출구 경계의 영향이 있는 전체 900 µm 적분과 인접 구간의 차이를 확인할 수 있습니다. ±500 µm 단면은 직선 채널의 전단응력 대조에 사용합니다. |
| `vortex_metrics.csv`, `separatrix.csv`, `upstream_vortex.png` | 단면을 통과하는 면 유량을 적분하여 유선함수를 재구성합니다. 코너에 연결된 재순환 영역의 길이, Lᵥ/Wᵤ, 내부 속도 및 주류와의 경계 유선을 기록합니다. 유체 쪽에서 와류가 분해되지 않으면 해당 길이와 속도는 빈 값입니다. |
| `mesh_comparison.csv` | 같은 물리적 구간의 접근 유량, 벽응력, 분포의 L¹ 차이입니다. 전체 구간과 양 끝을 2 µm씩 제외한 구간을 모두 기록합니다. |

논문 Fig. 8의 입자 없는 µPIV는 0.5 mL/h입니다. PEO 와류 내부 속도 약 0.2–2.0 mm/s와 Gly75 와류의 관측 분해능을 함께 비교합니다. Fig. 5(c)의 0.3 mL/h, 135분 사진에서는 목 가장자리의 돌출된 응집체를 식별할 수 있습니다. 이 투영 사진에서 최초 부착점이나 단차면의 점유 폭을 정량화하기 어려우므로 0.3 mL/h 그림에는 목 끝을 **정성적 기준점**으로 표시합니다. 추후 보정된 사진에서 측정한 구간은 `tools/compare_results.py`의 `--deposit-band LO HI`로 넣을 수 있습니다(µm).

## 전단점도를 일치시킨 탄성 민감도

두 번째 모드의 λ를 바꿀 때 ελ²를 고정하면 이 모델의 정상 단순 전단점도가 모든 전단률에서 보존됩니다. 각 모드의 지수 인자 f는 f² ln(f)=2ε(λγ̇)²를 만족하고 해당 모드의 점도 기여는 ηₚ/f입니다. 이 조건으로 첫 번째 모드와 CaBER에 대응하는 최장 시간은 유지하면서 탄성 민감도를 확인할 수 있습니다. 아래의 30 ms는 민감도 비교를 위한 예시값입니다.

```bash
python3 tools/shear_matched_fit.py --lambda2-ms 30 --write-config peo_lambda2_30ms.json
sbatch run_dep.sbatch peo 0.5 coarse --config peo_lambda2_30ms.json
```

새 피팅 파일은 별도 결과 경로를 사용합니다. 단독으로 `python3 tools/shear_matched_fit.py`를 실행하면 논문의 Cross 점도와 기본 ePTT의 정상 전단점도를 비교합니다.

## 초기화와 수렴 기준

ePTT 작업은 목표 유량의 10%에서 Cross 정상장을 구한 뒤, 같은 격자의 내부 속도·압력을 첫 점탄성 계산에 사용합니다. ePTT를 목표의 10%, 30%, 60%, 100%로 올리며, 각 단계가 수렴하면 다음 단계로 속도, 압력, 고분자 응력과 로그 conformation을 전달합니다. 각 단계의 입구는 그 유량에 맞춰 다시 생성합니다. 기본 응력 대류는 1차 풍상입니다.

기본값은 10회 반복마다 검사하며, U·p·두 모드의 정규화 초기 잔차가 10⁻⁵ 미만이어야 합니다. 유입·유출 차이와 지정 유량의 차이는 q의 10⁻⁵ 미만, 각 벽의 접근유량·속도분포·응력분포 변화는 10⁻³ 미만이어야 합니다. 출구의 최대 traction 오차도 출구 최대 여분응력의 10⁻³ 미만이어야 합니다. 3회 연속 통과하면 저장하고 종료합니다. 단계당 12,000회는 미수렴 계산의 상한입니다. 횟수와 계산 시간은 각 단계의 `run_status.json`, 후처리 결과에 쓰인 이산화와 물성은 `comparison.json`에서 확인할 수 있습니다.

격자는 `smoke` 807셀, `coarse` 6,552셀, `medium` 14,796셀, `fine` 33,510셀입니다. 단차와 목 주변 및 자유출구에서 세분화합니다. 출구 근처의 마지막 셀 길이는 기본 격자에서 약 0.73 µm입니다. `smoke`는 실행 경로 점검용으로, 수축부 와류가 더 촘촘한 격자와 다르게 분해될 수 있습니다.

## 이전 파일

이전 3D ParaView preset `tools/paraview_compare.py`는 2D 비교에서 사용하지 않으므로 삭제해도 됩니다. 같은 이름인 `run_dep.sbatch`와 `postprocess_dep.sbatch` 등은 overlay가 교체합니다. 이전 업데이트를 별도로 적용한 클러스터에만 있는 `run_newtonian_converged.sbatch`, `run_peo_diagnostic.sbatch`, `build_dep_update.sbatch` 및 그 보조 스크립트 `tools/configure_newtonian_stop.py`, `tools/collect_peo_diagnostics.py`는 이전 3D 작업의 재현용입니다. 새 실행 경로는 `runs/steady2d/`입니다.
