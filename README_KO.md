# Deposition: Gautschi / Foundation OpenFOAM-8

실험의 Gly75 Newtonian 유체와 PEO1.0의 2모드 exponential PTT 유체를 동일한 3D 수축 채널에서 비교하는 실행 패키지입니다. 입자를 해상하거나 부착·탈착을 계산하지 않습니다. 초기 입자 없는 유동장, 고분자를 포함한 벽 접선응력, 벽에 가까운 가상면의 대류 접근유량을 계산합니다.

## 1. 설치와 첫 실행

ZIP 내부의 최상위 폴더는 `Deposition/`입니다. `/home/lee6456`에서 압축을 풀면 원하는 위치가 됩니다.

```bash
cd /home/lee6456
unzip Deposition_OF8_rheoTool.zip
cd Deposition
sbatch build_dep.sbatch
```

빌드 로그 끝의 `BUILD COMPLETE`를 확인한 뒤:

```bash
sbatch run_dep.sbatch gly75 0.5 coarse
sbatch run_dep.sbatch peo 0.5 coarse
```

두 유체를 한 번에 제출하는 동등한 명령은 다음과 같습니다.

```bash
bash submit_comparison.sh 0.5 coarse
```

빌드 job ID가 예를 들어 `12345`라면, 빌드가 성공한 뒤에만 두 작업이 시작되도록 예약할 수도 있습니다. `12345`는 실제 제출 결과로 바꿉니다.

```bash
bash submit_comparison.sh 0.5 coarse --after-build 12345
```

실행 환경만 빠르게 확인하려면 `coarse` 대신 `smoke`를 사용합니다. smoke는 9,440셀이고 기본 0.0001초만 계산하며 두 번 저장합니다. **정상 유동의 물리적 비교에 쓰는 격자가 아닙니다.** 정상상태 조건을 충족하지 않으므로 `unconverged`가 예상되는 정상적인 검사 결과입니다.

## 2. 기존 OpenFOAM을 보존하는 방식

`config/environment.sh`의 기본값:

```bash
DEP_OF_ROOT=/home/lee6456/OpenFOAM/OpenFOAM-8
DEP_MODULES='gcc/11.4.1 openmpi/4.1.6'
```

이 module 조합은 이전 Gautschi OpenFOAM-8 작업에서 사용한 값입니다. 현재 설치를 다른 compiler/MPI 조합으로 빌드했다면 **그 조합으로 이 파일만 수정**합니다. 최신 compiler로 임의 교체하지 않습니다. module load 또는 기존 `blockMesh -help`가 실패하면 빌드는 중단합니다.

기본 Slurm account는 이전 Gautschi 작업의 `lilly-nan6`, partition은 `cpu`입니다. 다른 allocation을 쓴다면 파일 본문을 수정하거나 `sbatch -A YOUR_ACCOUNT ...`로 덮어쓸 수 있습니다. MPI 작업은 기본 16 ranks, rank당 1 CPU, 24시간입니다. 노드의 메모리 배정 정책에 맞춰 rank당 2 GB를 요청합니다.

다음은 모두 새 프로젝트 안에만 생성됩니다.

| 위치 | 내용 |
|---|---|
| `vendor/rheoTool/` | 고정 commit의 원본 구성방정식 소스와 라이선스 |
| `deps/` | SHA256로 검증하는 Eigen 3.4.0 원본 archive와 headers |
| `build/rheoToolSubset/` | OpenFOAM-8용으로 준비된 최소 라이브러리 소스 |
| `install/<WM_OPTIONS>/bin/` | `rheoDepositionFoam`, 물성 검증 실행 파일 |
| `install/<WM_OPTIONS>/lib/` | `libdepositionConstitutiveEquations.so` |
| `.venv/` | 결과 plot용 Python 환경 |
| `runs/` | 유체·유량·격자별 독립 계산 |
| `comparison/` | 비교 표, plot, ParaView state |

`OpenFOAM-8`, `ThirdParty-8`, `lee6456-8`, `$HOME/.bashrc`에 파일을 쓰거나 `Allwmake`를 실행하지 않습니다. OpenFOAM의 기존 라이브러리와 utilities만 참조합니다. 목적에 필요한 RheoTool 소스와 Eigen archive가 포함되어 있어 이 두 항목을 서버에서 다시 다운로드할 필요도 없습니다. Python plot 패키지 설치에는 네트워크가 필요합니다. 그 설치가 실패해도 CFD와 `.foam` 파일은 사용할 수 있습니다.

## 3. RheoTool 호환성과 구현 범위

공식 RheoTool은 OpenFOAM 7/9용 디렉터리를 제공하며, Foundation 8은 공식 테스트 대상에 포함되어 있지 않습니다. 따라서 이 패키지는 stock RheoTool 전체 설치가 아니라, 아래 commit의 필요한 부분을 OpenFOAM-8에서 실제 컴파일한 **프로젝트 전용 최소 포팅**입니다.

- 원본: https://github.com/fppimenta/rheoTool
- commit: `81f3e8cf6079e84afd04c5d847389c338404b268`
- constitutive source: `of90/src/libs/`
- static segregated solver: RheoTool의 SIMPLEC 알고리즘 기반
- executable: `rheoDepositionFoam`

구성방정식은 원본 `PTTLog`, `multiMode`, `Newtonian`, `CarreauYasuda`를 사용합니다. PETSc/Hypre의 block-coupled 선형 풀이, 이동격자, 전기유동, 열전달, 수동 스칼라는 이 실행 파일의 기능에 포함되지 않습니다. 해당 기능을 요청하면 지원하지 않는다고 명시적으로 중단합니다.

원본 대비 변경 사항은 `build/rheoToolSubset/PROVENANCE.json`에 남깁니다.

1. 사용하지 않는 block-coupled 의존성을 제거하고 OpenFOAM 기본 선형 솔버를 사용합니다.
2. log-conformation의 중복 고유값 한계(A=I 등)를 올바르게 처리합니다. 처음 응력이 0인 순수 전단에서도 theta_xy가 생성되어야 합니다. 이 수정은 구성방정식의 물성을 변경하는 것이 아니라 수치적 특이점의 올바른 한계 처리입니다.
3. 선택적 thermo dictionary가 없을 때의 null 접근을 제거합니다. 생성 케이스는 Constant thermo dictionary도 명시합니다.

원본 소스와 라이선스의 hash manifest를 검증하므로 검토되지 않은 upstream 업데이트를 자동으로 가져오지 않습니다.

## 4. 실험 조건과 유체 계수

`config/experiment.json`을 단일 입력 원본으로 사용합니다. 길이·유량·응력은 SI 단위입니다.

| 항목 | 값 |
|---|---|
| 상류 폭, 하류 폭 | 180 μm, 45 μm |
| 높이 | 100 μm |
| 상류 길이, 하류 길이 | 10 mm, 15 mm |
| 형상 | abrupt 4:1 contraction, sharp corner |
| 온도 | 20°C, 등온 가정 |
| 유량 | 0.3 / 0.5 / 0.7 mL/h |
| Gly75 | eta=0.035 Pa s, rho=1183 kg/m³ |
| PEO | rho=998 kg/m³ |

좌표는 x가 유동방향, 수축부 x=0, y가 폭, z가 높이입니다. 높이는 −50~50 μm입니다. 모든 벽은 no-slip, 입구는 지정 체적유량, 출구는 gauge pressure=0입니다. 입구 reservoir의 고분자 응력은 0으로 두고 실제 길이의 상류 채널 안에서 응력이 발달하도록 합니다. 수축부 바로 앞에 균일한 전단응력을 임의로 부여하지 않습니다.

두 모드의 ePTT:

| 모드 | eta_p (Pa s) | lambda (s) | epsilon | eta_s 기여 (Pa s) |
|---|---:|---:|---:|---:|
| 1 | 0.04112164164085602 | 0.04000000000000000 | 0.0008138270262392983 | 0.00695 |
| 2 | 0.15192835835914398 | 0.023488735159794105 | 0.8317237247424983 | 0 |

`zeta=0`, `destructionFunctionType exponential`입니다. RheoTool은 모드별 etaS를 더하므로 **0.00695를 두 모드에 중복 입력하지 않습니다.** 합계 eta0=0.2 Pa s입니다. 40 ms는 CaBER 값을 장시간 모드에 배정한 것이고 나머지 계수는 선택한 피팅입니다. 표의 많은 자릿수는 실험적 유효숫자가 아닙니다.

선택적인 `cross` 대조군도 있습니다. 실험의 세 번째 유체가 아니라 PEO의 탄성만 제외한 계산 대조군입니다. Carreau–Yasuda의 a=0.718, n=0.282를 사용하면 (n−1)/a=−1이므로 논문의 Cross 식과 정확히 같습니다.

```bash
sbatch run_dep.sbatch cross 0.5 coarse
```

## 5. 정상해를 찾는 방식과 계산량

이 실행 파일은 Euler 시간 적분과 SIMPLEC 내부 반복으로 유동·응력을 함께 풀어 정상상태를 찾습니다. 정상 방정식에서 없어지지 않는 응력 수송·변형 항을 유지합니다. 입자 궤적이나 부착 kinetics를 시간 적분하는 작업은 없습니다.

기본 maxCo=0.25, maxDeltaT=5e−5 s, 내부 반복 3회이며, 실제 timestep은 throat 격자와 속도에 따라 더 작아집니다. 기본 시간 상한 4초는 계산 완료의 물리적 보증이 아닙니다. 고 Wi에서 비정상 유동이 유지되면 정상상태로 표시하지 않습니다.

정상상태 판정은 최소 0.4초 이후 다음 조건을 0.2초 동안 유지하는 것입니다. 기본 0.01초 간격의 진단과 저장 시점에서 검사합니다.

- 수축부 주변 |x|<180 μm에서 U와 전체 extra stress의 체적 가중 L2 상대 변화 <1e−4.
- step-face 평균 접선응력, 대류 접근유량, 압력강하의 상대 변화 <1e−4.
- 입구·출구 유량 불균형 <1e−5.

조건을 충족하면 최종 field와 `stationary.json`을 저장하고 자동 종료합니다. 조건을 충족하지 않은 시간 상한 도달은 `unconverged`입니다. `run_status.json`에서 `stationary/unconverged/failed/interrupted`를 구분합니다. 선형 residual이 작거나 Slurm job이 완료됐다는 것만으로 정상해라고 판단하지 않습니다.

| 격자 | 셀 수 | 용도 |
|---|---:|---|
| smoke | 9,440 | 설치·입출력 검사, 물리적 결과 해석 금지 |
| coarse | 268,992 | 첫 정상해·동작 확인 |
| medium | 925,488 | 정량 비교의 다음 단계 |
| fine | 3,119,472 | 격자 독립성 검사 |

처음에는 0.5 mL/h coarse의 두 유체를 비교하고, 같은 유량의 medium 결과와 벽응력·접근유량을 비교합니다. sharp corner의 최대 응력 하나를 격자 수렴 대상으로 삼지 않습니다. 시간 적분 오차도 maxCo를 절반으로 줄여 확인해야 합니다. coarse 계산이 빠르게 끝난다고 보장하는 패키지는 아니며, 실제 steady까지의 소요시간은 Gautschi에서 확인해야 합니다.

## 6. 출력과 정확한 물리량 정의

예를 들어 PEO 출력:

```text
runs/peo_Q0p5_coarse/
  case.foam
  metadata.json
  run_status.json
  postProcessing/deposition/
    metrics.csv
    approach_latest.csv
    stationary.json           # 조건 충족 시에만 존재
```

주요 field:

| 이름 | 단위 및 의미 |
|---|---|
| `U` | 속도, m/s |
| `p` | kinematic pressure, m²/s² |
| `p_Pa` | rho*p, Pa |
| `tau` | 두 모드 고분자 응력의 합, Pa; Newtonian에는 없을 수 있음 |
| `tauTotal` | 고분자+배경점성 extra stress, Pa; 압력 제외 |
| `N1` | 좌표축 기준 polymer tau_xx−tau_yy, Pa |
| `wallShearStressTotal` | 전체 extra stress의 벽 접선 traction 벡터, Pa |
| `theta_mode*` | log conformation, 무차원 |
| `eigVals_mode*`, `eigVecs_mode*` | 재시작에 필요한 conformation 고유값·고유벡터 |

일반 수축유동에서 `N1`은 지정 좌표축의 정상응력 차이이며, 국소 흐름에 맞춰 회전한 단순전단 물성 N1과 항상 같지는 않습니다. 각 field는 실제 계산값이며 전단 물성곡선의 값을 위치마다 대신 넣지 않습니다.

`wallShearStressTotal = (I−nn) · tauTotal · n`입니다. 벽의 외향 법선에 대한 traction을 출력하며 비교 plot은 그 크기를 사용합니다. 압력은 접선 투영에서 제거됩니다. 일반 OpenFOAM의 점성 부분만 포함한 wallShearStress 값을 사용하는 것이 아닙니다.

`wallShearMean_Pa`는 **상류 step face**의 면적 가중 평균입니다. 안쪽 throat edge, 바깥쪽 edge, 높이방향 edge에서 각각 2 μm를 제외합니다. `downstreamWallShearMean_Pa`는 하류 측벽의 x=2~180 μm 구간 평균으로 같은 높이 edge 제외를 사용합니다. 최대값은 참고값입니다.

`approachFlux_m3_s`는 step face에서 2 μm 떨어진 x=−2 μm 면의 양의 Ux 적분입니다. 두 wing에 대해 24×32개의 midpoint sampling을 사용하고 OpenFOAM의 cell-point interpolation으로 속도를 구합니다. 적분 범위는 24.5<|y|<88 μm, |z|<48 μm입니다. `approachNetFlux_m3_s`는 양·음 방향을 상쇄한 순유량도 함께 제공합니다. `approach_latest.csv`에 모든 sample의 위치·속도·면적이 기록됩니다.

이 값은 **유체의 대류 접근량**으로, 불투과 벽 자체의 유입유량이나 실제 포집률이 아닙니다. 균일한 입자 수농도 c0를 별도로 정하면 c0×approachFlux가 입자수 접근량의 추정치가 됩니다. 이는 유선을 반복 통과하는 입자를 중복 집계할 수 있는 gross flux이며 실제 부착률은 아닙니다. 공간격자와 sampling resolution을 모두 확인해야 합니다. `constant/depositionProperties`의 sampleNy/sampleNz를 두 배로 바꿔 적분의 수렴을 확인할 수 있습니다.

같은 물리적 거리의 값을 비교해야 하므로 offset을 첫 번째 cell centre 위치로 자동 변경하지 않습니다. 다른 offset은 새 case의 `--approach-offset-um` 입력으로 만들 수 있습니다. 저수준 case 생성은 `python3 tools/make_case.py --help`를 참조합니다.

## 7. 비교 plot과 ParaView

두 작업이 끝난 후:

```bash
sbatch postprocess_dep.sbatch 0.5 coarse
```

`comparison/Q0p5_coarse/`에 비교 CSV/JSON과 응력·유량·압력 이력 PNG를 만듭니다. 서로 다른 시점에서 정상상태에 도달해도 각 유체의 최종 결과를 비교하며, 두 marker가 모두 일치할 때만 `STEADY`로 표시합니다. 그 외에는 `UNCONVERGED PREVIEW`입니다.

서버에서 pvpython이 사용 가능하면 ParaView state도 만듭니다. headless 서버에서는 기본적으로 state만 생성하므로 X server가 필요하지 않습니다. 설치된 pvpython 경로는 `DEP_PVPYTHON`으로 지정할 수 있습니다. 특정 Gautschi ParaView module 이름을 임의로 가정하지 않습니다.

로컬 PC의 ParaView에서는 계산 폴더의 `case.foam`을 열 수 있습니다. 같은 PC에서 ParaView의 pvpython으로 아래를 실행하면 공통 색상 범위의 velocity/streamline, N1, wall traction PNG와 편집 가능한 `.pvsm`을 만듭니다.

```bash
pvpython tools/paraview_compare.py \
  runs/gly75_Q0p5_coarse runs/peo_Q0p5_coarse \
  --output comparison/Q0p5_coarse/paraview
```

화면이 없는 환경에서는 끝에 `--state-only`를 붙입니다. Windows에서는 설치된 ParaView의 `bin/pvpython.exe` 경로로 실행합니다. 원격에서 PC로 옮길 때 reconstructed 최신 time, `constant`, `system`, `metadata.json`, `case.foam`, `postProcessing/deposition`을 같이 가져오면 processor 폴더 없이 볼 수 있습니다. 재시작하려면 서버의 processor 폴더는 보존합니다.

## 8. 재시작과 다른 유량

기존 case는 자동 삭제하거나 다시 생성하지 않습니다. 다음 명령은 저장된 완전한 checkpoint에서 이어서 계산합니다.

```bash
sbatch run_dep.sbatch peo 0.5 coarse --resume --end-time 8
```

MPI rank 수는 기존 case와 같아야 합니다. 모든 rank에서 U/p/tau/theta/eigenstate가 완전히 기록된 최신 공통 시간을 사용합니다. 이전 정상상태 marker와 rollback 시점 뒤 history는 별도 이름으로 보존합니다. Slurm 제한시간 120초 전에는 현재 상태 저장을 요청합니다. 급작스러운 노드 장애 시 마지막 완전한 주기 checkpoint를 사용합니다.

다른 유량의 독립 명령:

```bash
sbatch run_dep.sbatch gly75 0.3 coarse
sbatch run_dep.sbatch peo 0.3 coarse
sbatch run_dep.sbatch gly75 0.7 coarse
sbatch run_dep.sbatch peo 0.7 coarse
```

정상 실행 중 시간 상한·write 간격 등만 변경하려면 runtime dictionary를 수정할 수 있지만, constitutive parameters·유량·격자를 바꾸어 기존 결과를 같은 케이스에 섞지 않습니다. 모델 수정 후에는 새 case를 생성합니다.

## 9. 검증 범위와 원고와의 비교

이 패키지의 로컬 검증 결과는 `VALIDATION.md`와 `validation/`의 작은 evidence 파일을 참조합니다. Gautschi 계정에 직접 접속하거나 작업을 제출한 것은 아닙니다. 최초 cluster build는 실제 module/ABI 환경에서의 마지막 확인 단계입니다. source/toolchain/library fingerprint가 달라지면 run script가 재빌드를 요구합니다.

원고의 기존 CFD는 **2D, 상·하류 각각 0.9 mm**의 영역이었습니다. 여기서는 사용자가 요청한 실험 형상의 **3D, 10 mm+15 mm**를 사용하므로 원고의 2D 유동장과 완전히 같은 계산은 아닙니다.

또한 원고에 기재된 단면·길이·Gly75 점도를 그대로 사용한 rigid no-slip 직사각 duct의 완전발달 압력강하는 수축손실을 제외해도 Q=0.3/0.5/0.7에서 약83.3/138.9/194.5 kPa입니다. 원고 Table1의38.4/63.5/84.2 kPa와 차이가 있어, Table1 압력 일치 여부를 이 코드의 통과 기준으로 사용하지 않습니다. 형상이나 점도를 이 차이에 맞추려고 임의 조정하지 않았습니다.

공식 자료:

- RheoTool: https://github.com/fppimenta/rheoTool
- Foundation OpenFOAM-8: https://openfoam.org/version/8/
- Gautschi queues: https://docs.rcac.purdue.edu/userguides/gautschi/run_jobs/queues/
- Gautschi module/toolchain overview: https://docs.rcac.purdue.edu/userguides/gautschi/using_ai_agents/

GPL 원본/수정 소스의 배포 조건은 `COPYING` 및 각 파일의 license header를 따릅니다.
