# P2RANK Installation Guide for Windows

## Step 1: Install Java

P2RANK requires Java 8 or later.

### Option A: Using Chocolatey (Recommended)
```bash
# Install Chocolatey package manager if not installed
# https://chocolatey.org/install

# Then install Java
choco install openjdk11
```

### Option B: Manual Download
1. Download OpenJDK from: https://adoptium.net/
2. Choose: Latest LTS (Java 17 or 21)
3. Install and add to PATH

### Option C: Using Scoop
```bash
scoop install openjdk
```

## Step 2: Verify Java Installation
```bash
java -version
```

Should output something like:
```
openjdk version "11.0.x" or higher
```

## Step 3: Download P2RANK

```bash
cd C:/Users/Jahyun/PycharmProjects/ASMC/pdb_structures

# Download P2RANK (latest version 2.4.2)
curl -L https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz -o p2rank.tar.gz

# Extract
tar -xzf p2rank.tar.gz

# or using 7-Zip if tar doesn't work:
# 7z x p2rank.tar.gz
# 7z x p2rank.tar
```

## Step 4: Run P2RANK

```bash
cd pdb_structures

# Windows - use prank.bat
./p2rank_2.4.2/prank.bat predict pdb3rfv_chainA.pdb -o p2rank_results

# Or if the above doesn't work:
java -jar p2rank_2.4.2/p2rank.jar predict pdb3rfv_chainA.pdb -o p2rank_results
```

## Step 5: Compare Results with Known Active Sites

```bash
# Run our comparison script
python compare_p2rank_with_known.py p2rank_results/pdb3rfv_chainA.pdb_predictions.csv
```

## Expected P2RANK Performance

Based on P2RANK's benchmark performance:
- **Expected F1-score: 50-70%** (vs our 7-19%)
- Much better at identifying functional sites
- Uses machine learning trained on 15,000+ structures

## Alternative: Use Online P2RANK

If local installation fails, use the web server:
https://prankweb.cz/

1. Upload `pdb3rfv_chainA.pdb`
2. Download prediction results
3. Compare with our known active sites
