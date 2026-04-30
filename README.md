# 🧬 Ann2Bool: scRNA-seq to Boolean Network Pipeline

단일 세포 전사체 데이터(scRNA-seq)를 기반으로 전사 인자(TF) 중심의 유전자 조절 네트워크(GRN)를 추론하고, 이를 Boolean Network 모델로 구축하고 타겟을 발굴하는 통합 프레임워크입니다. 정상 세포에서 특정 암(예: PDAC) 등으로 진행되는 궤적 상의 핵심 조절 인자를 발굴하고, in-silico Perturbation 분석을 통해 세포 운명 제어 타겟을 탐색할 수 있습니다.

## 🌟 Key Features
* **Data Integration:** 단일 세포 발현량(AnnData)과 Monocle3 기반의 Pseudotime 궤적 데이터 통합
* **GRN Inference:** `pySCENIC`을 활용한 세포 타입 특이적 TF-TF 상호작용 및 조절 방향(Activation/Inhibition) 추론
* **Automated Rule Learning:** GMM(Gaussian Mixture Model)을 통한 TF 활성도 이진화 및 시간 흐름 기반 Boolean 논리 규칙 자동 학습 (Decision Tree / QMC 적용)
* **Perturbation Analysis:** 구축된 `.bnet` 모델에서 구해지는 Attractor를 실제 세포 상태에 매핑하고, TF Knockout/Overexpression 시 특정 세포 타입으로의 전이 확률 변화 분석

## 🛠 Prerequisites

대규모 전사 인자 스크리닝 및 네트워크 추론(GRNBoost2 등) 과정에서 많은 컴퓨팅 리소스가 요구되므로, 충분한 메모리가 확보된 **Linux 서버 및 HPC 환경**에서의 실행을 권장합니다. 

이 파이프라인은 **Python 3.10.20**을 기반으로 작성되었습니다. 환경 충돌을 방지하기 위해 제공되는 `ann2bool_environment.yml` 파일을 사용하여 동일한 Conda 가상환경(`ann2bool`)을 간편하게 구축할 수 있습니다.

**Conda 가상환경 생성 및 패키지 설치:**
```bash
# 1. 터미널에서 리포지토리를 다운로드(Clone)한 폴더로 이동합니다.
# cd path/to/Ann2Bool

# 2. 제공된 yml 파일을 사용하여 'ann2bool' 가상환경을 생성하고 패키지를 설치합니다.
conda env create -f ann2bool_environment.yml

# 3. 설치가 완료되면 생성된 가상환경을 활성화합니다.
conda activate ann2bool
```
## 📥 필수 입력 데이터 (Required Inputs)
본 파이프라인은 다음과 같은 형식의 데이터를 인풋으로 사용합니다:

* **1. AnnData 객체 (.h5ad)**: 
    - HVG (Highly Variable Genes) 로 subset이 완료된 상태여야 합니다.

    - Raw Counts: `adata.layers['counts']` 와 같이 원본 count 데이터가 포함되어야 합니다.

    - Cell Type: `adata.obs` 내에 세포 타입 정보가 포함되어야 합니다.

* **2. Pseudotime data**: 각 세포별 Pseudotime 값이 저장된 별도의 .csv 파일이 필요합니다.
    - 첫 번째 열에 세포 ID(cell_id), 두 번째 열에 Pseudotime 값(PT)을 갖는 형식이어야 합니다.

        | cell_id | PT |
        |---|---|
        | cell1 | PT value 1 | 
        | cell2 | PT value 2 |
        | cell3 | PT value 3 |
        | ... | ... |


## 🛠️ 주요 기능 및 워크플로우

| 단계 | 명칭 | 설명 |
| :--- | :--- | :--- |
| **STEP 1** | **데이터 전처리** | AnnData 로드, Monocle3 의사시간 병합, 타겟 세포 타입 선택 및 정규화 수행|
| **STEP 2** | **pySCENIC GRN 추론** | GRNBoost2 및 cisTarget을 사용하여 Regulon을 추론하고, RSS(Regulon Specificity Score)로 핵심 TF 선별 |
| **STEP 3** | **TF-TF 네트워크 추출** | Regulon 내 TF 간 AUC 상관관계를 분석하여 활성(+) 및 억제(-) 엣지(Edge) 판별 |
| **STEP 4** | **이진화 (Binarization)** | GMM(Gaussian Mixture Model)을 사용하여 TF 활성 스코어를 0(Off)과 1(On) 상태로 변환 |
| **STEP 5** | **Boolean rule 학습** | **Quine-McCluskey(QMC)** 또는 **Decision Tree** 방식을 선택하여 $t \to t+1$ 전이 규칙 도출 |
| **STEP 6** | **네트워크 정제 및 내보내기** | 고립된 상수 노드 및 무의미한 Self-loop를 제거하고 `.bnet` 파일로 저장 |
| **STEP 7** | **attractor-cell 매핑** | 계산된 어트랙터를 Euclidean distance 기반으로 실제 UMAP 공간의 세포들과 매핑 |
| **STEP 8** | **Perturbation 및 시각화** | 유전자 Knock-out/Overexpression 시뮬레이션 결과 분석 및 네트워크 토폴로지 시각화 |

## ⚙️ 설정 (Configuration)
`CONFIG` dictionary 부분을 통해 분석 파라미터를 조정할 수 있습니다:
*   `rule_method`: 불리언 규칙 학습 알고리즘 선택 (`"qmc"` 또는 `"decision_tree"`)
*   `rss_threshold`: 세포 타입별 특이적 TF 선별 기준 (기본값: 0.3)
*   `window_size` / `window_step`: QMC 학습 시 시계열 데이터 노이즈 완화를 위한 슬라이딩 윈도우 설정

## 📂 주요 출력물
*   `pdac_boolean_network.bnet`: 최종 도출된 불리언 네트워크 모델 파일
*   `auc_matrix.csv`: 세포별 TF 활성도 행렬
*   `cellstate_binary_profile.csv`: 세포 타입별 대표 이진 상태 프로필
*   `attractor_umap.pdf`: UMAP 상에 매핑된 어트랙터 위치 및 Basin 크기 시각화
*   `perturbation_basin_summary.csv`: 각 유전자 섭동에 따른 세포 운명 변화 요약

---

**Note**: 본 파이프라인은 특히 ADM(Acinar-to-Ductal Metaplasia)에서 PDAC(Pancreatic Ductal Adenocarcinoma)로의 전환과 같은 암 발생 기전 분석을 기반으로 설계되었습니다.
