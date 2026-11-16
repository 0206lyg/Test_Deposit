# Newtonian 수렴 종료 + PEO 오류 진단 업데이트

기존 OpenFOAM-8 설치와 실험 물성은 변경하지 않습니다. 이 업데이트는 기존 Deposition 패키지에 적용하는 소스/실행 스크립트 묶음입니다. Newtonian 종료 방식을 바꾸고, PEO의 실패 원인을 기록합니다. PEO 불안정성을 고쳤다고 주장하는 업데이트가 아닙니다.

## 1. 기존 Newtonian 작업이 아직 실행 중이면 저장 후 종료

실행 중인 solver가 사용하는 라이브러리를 덮어쓰지 않도록 **기존 작업이 끝난 뒤 빌드**합니다. 현재 Newtonian 로그에서는 유동장과 관측량이 이미 안정되어 있으므로, 다음 명령으로 현재 상태를 저장하고 정상적으로 solver를 종료할 수 있습니다.

```bash
cd /home/lee6456/Deposition
source tools/load_environment.sh
foamDictionary /scratch/gautschi/lee6456/results_deposition/gly75_Q0p5_coarse/system/controlDict -entry stopAt -set writeNow
```

`log.solver.*`의 `End`와 Slurm 작업 종료를 확인합니다. 기존 실행 스크립트가 결과를 복원할 시간을 줍니다. 이전 버전에서는 이 수동 종료가 `unconverged`로 기록될 수 있습니다. 새 버전으로 아래 재시작을 하면 새로운 수렴 조건으로 다시 확인하고 저장합니다. 이미 종료된 작업에는 위 명령이 필요 없습니다.

## 2. 업데이트 설치 및 빌드

ZIP을 `/home/lee6456/Deposition`으로 옮긴 뒤:

```bash
cd /home/lee6456/Deposition
unzip Deposition_update_20261116.zip
python3 Deposition_update_20261116/apply_update.py /home/lee6456/Deposition
sbatch build_dep_update.sbatch
```

`UPDATE BUILD COMPLETE`가 나온 뒤 계산을 제출합니다. 기존 `run_dep.sbatch`, `build_dep.sbatch`, `config/environment.sh`와 케이스는 보존합니다. 따라서 사용자가 바꾼 계정/코어 수/경로도 그대로 남습니다. 덮어쓰는 기존 소스는 `update_backups/<timestamp>/`에 백업합니다. 예상과 다른 소스 수정이 있으면 아무것도 덮어쓰지 않고 멈춥니다.

새 sbatch 세 개는 첨부된 실행 스크립트 기준 계정 `lilly-roa3`, partition `cpu`를 사용합니다. 실행은 4 MPI ranks, 빌드는 8 CPU입니다. 필요한 경우 각 파일의 SBATCH 줄만 변경합니다.

## 3. Newtonian: 기존 결과에서 재시작 (권장)

기존 작업이 종료되었고 체크포인트가 있다면:

```bash
sbatch run_newtonian_converged.sbatch 0.5 coarse --resume --case /scratch/gautschi/lee6456/results_deposition/gly75_Q0p5_coarse
```

이 명령은 기존 케이스의 종료 설정을 백업한 뒤 새 조건으로 바꾸고, 가장 최근의 완전한 4-rank 체크포인트에서 재시작합니다. 격자/물성/유량/코어 수를 바꾸지 않습니다. 재시작 후 수렴 비교 기준은 새로 수집합니다.

새 계산으로 검증하려면:

```bash
sbatch run_newtonian_converged.sbatch 0.5 coarse
```

새 폴더는 `/scratch/gautschi/lee6456/results_deposition/gly75_Q0p5_coarse_converged`입니다. 이미 존재하면 덮어쓰지 않고 종료합니다.

### 변경된 Newtonian 종료 조건

- `minimumSteadyTime=0`, `steadyWindow=0`: 0.4초/0.2초를 강제로 기다리지 않습니다.
- 진단 간격은 물리 시간 0.01초이며, 서로 독립적인 비교를 **3회 연속** 통과하면 즉시 저장·종료합니다.
- 첫 진단은 비교 기준을 저장하며 통과 횟수로 세지 않습니다.
- 수축부 ROI의 속도·전체 응력 RMS 변화, 수축면 평균 벽 접선응력, 가상면 접근유량, 압력강하, 질량보존을 함께 확인합니다.
- 변화 조건: `절대 변화 <= 절대 허용오차 + 상대 허용오차 × 크기`. 상대 오차는 1e-4, 질량보존은 1e-5입니다. 절대 허용오차는 속도 1e-9 m/s, 응력/벽응력 1e-6 Pa, 압력 1e-3 Pa, 접근유량/유량 불일치 1e-20 m³/s입니다.
- 저장 시점의 추가 진단은 수렴 횟수를 늘리지 않습니다. 로그에 `stable checks=3/3`와 종료 사유를 기록합니다.
- `stationary.json`과 최종 fields, `run_status.json`이 일치하는지 기존 검사기로 확인합니다.

4초 endTime은 미수렴 계산의 상한으로 남으며, 수렴 시 그 시간까지 계산하지 않습니다. PEO의 정상상태 판정 정책은 이 업데이트에서 완화하지 않았습니다.

## 4. PEO: 원래 조건을 유지한 짧은 진단 실행

```bash
sbatch run_peo_diagnostic.sbatch 0.5 coarse
```

새 폴더: `/scratch/gautschi/lee6456/results_deposition/peo_Q0p5_coarse_diagnostic`

이전 실패 케이스를 덮어쓰지 않습니다. 원래 2모드 ePTT 계수, 격자, 최대 Courant 수, 시간 간격 제어, 내부 반복 수, 수치 scheme을 유지합니다. 계산 상한만 0.003초로 짧게 합니다. 원래 0.1초 저장 간격을 유지하여 출력 시각에 맞춘 시간 간격 조절이 재현 경로를 바꾸지 않도록 했습니다. 실패 시점의 field snapshot은 별도로 저장합니다. 원래 오류 시점은 약 0.001478초입니다. 0.003초를 통과하더라도 모든 장시간 조건이 안정적이라는 의미는 아닙니다.

진단 코드는 다음을 구분합니다.

1. 구성방정식 입력의 비유한 값 또는 양수가 아닌 conformation 고유값.
2. ePTT 완화계수 `f=exp[epsilon/(1-zeta)*(tr(C)-3)]`의 지수/곱셈 overflow.
3. theta 업데이트 뒤 conformation을 복원하는 지수함수의 문제.
4. conformation에서 응력으로 변환할 때의 overflow.

원래 코드가 실패하기 직전에 가능한 경우 모든 rank가 진단 snapshot을 저장하고 명시적으로 중단합니다. 어느 모드, 시각, 내부 반복, processor, 셀/경계면, 좌표에서 발생했는지와 theta/tau/C 고유값/속도구배/지수 인수를 기록합니다. 직전 시간의 theta와 직전 성공한 구성방정식 계산의 속도구배도 함께 보존합니다. 지수 인수 및 `log(deltaT*f/lambda)` 이력을 기록하여 최종 overflow 이전의 급격한 증가를 조사할 수 있게 합니다. 이 비율 자체가 완전한 비선형 안정성 기준인 것은 아닙니다.

**응력이나 지수 인수를 잘라내지 않고, epsilon/lambda를 바꾸지 않습니다.** 추가 검사는 계산비용이 있으므로 이 실행은 오류 진단용입니다. 로그의 모든 선형 residual이 작더라도 비선형 응력 업데이트의 안정성을 보장하지 않습니다.

## 5. 다음에 전달할 파일은 하나

성공 또는 실패 후 자동으로 생성합니다:

```text
/scratch/gautschi/lee6456/results_deposition/peo_Q0p5_coarse_diagnostic/diagnostic_bundle.tar.gz
```

설정 파일, 로그 앞/끝, `EPTT_DIAGNOSTIC` 이력, 실패 지점의 수치, 체크포인트 상태가 들어 있습니다. 전체 field snapshot은 큰 파일이므로 이 묶음에서는 제외하고 케이스의 `processor*/postProcessing/epttDiagnostics/`에 보관합니다. snapshot은 재시작용 디렉터리가 아니며 체크포인트로 선택되지 않습니다.

예상치 못한 강제 종료로 자동 수집이 실행되지 않았다면:

```bash
cd /home/lee6456/Deposition
python3 tools/collect_peo_diagnostics.py /scratch/gautschi/lee6456/results_deposition/peo_Q0p5_coarse_diagnostic
```

추가로 현재의 긴 로그 전체를 다시 전달할 필요는 없습니다. 위 진단 묶음으로 다음 원인 분석을 진행합니다.

## 검증 범위

동봉한 `VALIDATION_UPDATE.md`에 이번 버전에서 실제 수행한 빌드·실행·오류 검출 시험을 구분해 기록했습니다. Gautschi의 실제 오류 셀과 최초 원인은 새 진단 실행에서 확인해야 합니다.
