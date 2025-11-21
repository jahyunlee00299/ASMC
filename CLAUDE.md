# ASMC Project Context

## Project Overview
ASMC (Active Site Mining and Clustering) - 단백질 active site 기반 서열 클러스터링 도구

## Key Reference Structure: PDB 3RFV (AtUDH)

### PDB Numbering Offset
**PDB 3RFV는 residue 2번부터 시작합니다.**

```
PDB_num = alignment_pos + 2
alignment_pos = PDB_num - 2
```

**핵심 촉매 잔기:**
- **Y136 (PDB numbering)** = Y134 (alignment position)
- NAI (NADH)와 2.54Å, 15L (galacturonate)와 2.53Å 거리

### AlphaFold 구조와의 Offset 차이
AlphaFold 구조는 residue 1번부터 시작:
- 3RFV Y136 → AlphaFold ~Y134-140 (구조마다 다름)
- **서열 정렬 vs 구조적 거리**로 대응 잔기가 다를 수 있음
- 예: A0A0N8HGG8에서 서열정렬=Y139, 구조거리=Y140

## Active Site Definition Files

| 파일 | 설명 | Residues |
|------|------|----------|
| `atudh_active_site.txt` | NAD+ binding (Direct+Proximal) | 31 |
| `atudh_nad_site.txt` | NAD+ binding only | 31 |
| `atudh_substrate_site.txt` | Substrate binding (15L) | 12 |
| `atudh_combined_3_5.5A.txt` | Combined (3.0-5.5Å cutoff) | 40 |

위치: `data/pdb_structures/`

## Key Files

- `asmc/active_site_extractor.py` - ActiveSiteExtractor 클래스
- `test_data/AtUdh_pdb3rfv_chainA.pdb` - Reference PDB
- `test_data/UDHs_filtered_std2.5.fasta` - UDH 서열 데이터셋
- `docs/ActiveSiteExtractor.md` - 상세 문서

## Important Notes

1. Active site 추출 시 반드시 offset 확인
2. PyMOL에서는 PDB 번호 사용 (Y136)
3. 서열 정렬 기반 분석에서는 alignment position 사용 (Y134)
4. AlphaFold 구조는 별도 offset 고려 필요
