# Python 3.14 호환성 정보

## ✅ 호환성 상태
이 프로젝트는 **Python 3.14**에서 정상적으로 작동하도록 업데이트되었습니다.

## 📋 업데이트된 패키지 버전

### 기존 버전 (Python 3.8-3.13용)
```
biopython==1.85
numpy==2.2.5
scikit-learn==1.6.0
scipy==1.15.2
pillow==11.1.0
plotnineseqsuite==1.1.3
pytest==8.3.4
xlsxwriter==3.1.1
```

### Python 3.14 호환 버전 (업데이트 완료)
```
biopython==1.86
numpy==2.3.4
scikit-learn==1.7.2
scipy==1.16.3
pillow==12.0.0
plotnineseqsuite==1.2.0
pytest==8.4.2
xlsxwriter==3.2.9
```

## 🔧 수정된 파일들

1. **requirements.txt** - Python 3.14 호환 버전으로 업데이트
2. **pyproject.toml** - 의존성 최소 버전 업데이트
3. **asmc/__init__.py** - 모든 함수와 예외 클래스 export 추가
4. **tests/test_asmc.py** - Windows 경로 이스케이프 문제 수정
5. **tests/test_utils.py** - Windows 경로 이스케이프 문제 수정

## ✨ 주요 변경사항

### 1. NumPy 버전 업그레이드
- NumPy 2.2.5는 Python 3.14에서 컴파일 실패
- NumPy 2.3.4로 업그레이드하여 해결

### 2. 테스트 코드 수정
- Windows 파일 경로에서 백슬래시 문제 해결
- `re.escape()` 함수를 사용하여 정규식 패턴 이스케이프 처리

### 3. 모듈 Import 문제 해결
- `asmc.__init__.py`에 모든 함수 export 추가
- `RenumberResiduesError` 예외 클래스 export 추가

## 📊 테스트 결과
```
============================= test session starts =============================
platform win32 -- Python 3.14.0, pytest-8.4.2, pluggy-1.6.0
collected 30 items
tests/test_asmc.py ................... [19 passed]
tests/test_utils.py ............ [11 passed]
============================= 30 passed in 2.17s ==============================
```

## 🚀 설치 방법

### 개발 모드 설치 (권장)
```bash
pip install -e .
```

### requirements.txt 사용
```bash
pip install -r requirements.txt
```

## ⚠️ 주의사항

1. **Python 버전**: Python 3.14.0 이상 필요
2. **Windows 사용자**: Visual Studio Build Tools 필요 (NumPy 컴파일용)
3. **메모리**: 대용량 데이터셋 처리시 충분한 RAM 필요

## 🔄 이전 Python 버전으로 돌아가기
Python 3.13 이하 버전을 사용하는 경우, Git에서 이전 버전의 requirements.txt를 체크아웃하세요:
```bash
git checkout 9434035 -- requirements.txt pyproject.toml
```

## 📝 변경 이력
- **2025-11-04**: Python 3.14 호환성 업데이트 완료
- 모든 의존성 패키지를 최신 호환 버전으로 업데이트
- 30개 테스트 케이스 모두 통과

---
작성자: ASMC 유지보수팀
날짜: 2025년 11월 4일