# Uronate Dehydrogenase AlphaFold Structures

이 디렉토리는 AlphaFold로 예측된 UDH 구조들을 저장합니다.

## 필요한 파일

다음 5개 구조의 PDB 파일을 이 디렉토리에 복사하세요:

| UniProt ID  | Rossmann RMSD | 설명      | 파일명 예시 |
|-------------|---------------|-----------|------------|
| A0A1I2LZE5  | 15.78 Å       | 가장 유사  | AF-A0A1I2LZE5-F1-model_v4.pdb |
| A0A1I1AQL9  | 17.27 Å       | -         | AF-A0A1I1AQL9-F1-model_v4.pdb |
| A0A2Z4AB20  | 19.85 Å       | -         | AF-A0A2Z4AB20-F1-model_v4.pdb |
| A0A1N6JJV2  | 20.89 Å       | -         | AF-A0A1N6JJV2-F1-model_v4.pdb |
| A0A1U7CQA8  | 25.19 Å       | 가장 다름  | AF-A0A1U7CQA8-F1-model_v4.pdb |

## AlphaFold DB에서 다운로드

각 구조는 AlphaFold Protein Structure Database에서 다운로드할 수 있습니다:

```bash
# 예시: A0A1I2LZE5 다운로드
wget https://alphafold.ebi.ac.uk/files/AF-A0A1I2LZE5-F1-model_v4.pdb

# 또는 웹 브라우저에서:
# https://alphafold.ebi.ac.uk/entry/A0A1I2LZE5
```

## 파일명 요구사항

- 파일명에 **UniProt ID가 포함**되어야 합니다
- 예: `A0A1I2LZE5_model.pdb`, `AF-A0A1I2LZE5-F1-model_v4.pdb` 모두 OK
- 스크립트가 파일명에서 ID를 자동으로 찾습니다

## Windows에서 파일 복사

```cmd
REM 원본 폴더에서 이 폴더로 복사
copy C:\path\to\alphafold\structures\*.pdb .
```

## Linux/Mac에서 파일 복사

```bash
# 원본 폴더에서 이 폴더로 복사
cp /path/to/alphafold/structures/*.pdb .
```

## 확인

```bash
# 파일이 제대로 복사되었는지 확인
ls -lh *.pdb

# 각 파일명에 ID가 포함되어 있는지 확인
ls *.pdb | grep -E "A0A1I2LZE5|A0A1I1AQL9|A0A2Z4AB20|A0A1N6JJV2|A0A1U7CQA8"
```

모든 파일이 준비되면, 상위 디렉토리에서 `extract_and_align_active_sites.py`를 실행하세요!
