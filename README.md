# dna-sequence-analyzer

## Features
- Local alignment (Smith-Waterman algorithm)
- Global alignment (Needleman-Wunsch algorithm)
- Edit distance and similarity (Needleman-Wunsch algorithm)
- RNA to amino acids translation

## Available commands
```
Usage: analyze.py [OPTIONS] SEQUENCE_A SEQUENCE_B

Options:
  -S, --summary
  -s, --similarity
  -e, --edit-distance
  -a, --alignment [global|local]
  --load-csv                      Load scores.csv and edit_cost.csv
  --help                          Show this message and exit.
```
```
Usage: translate.py [OPTIONS] [SEQUENCE]

Options:
  -i, --input-file FILE  Path to text file containing long nucleotide sequences (1 sequence = 1 line)
  --help                 Show this message and exit.
PS D:\Repos\dna-sequence-analyzer>
```
## Requirements
- Python 3.7 (type annotations)
- numpy (storing matrices)
- pandas (loading CSV into DataFrame)
- click (CLI interface)

`pip install -r requirements.txt`

## Usage examples
```
python analyze.py AGCT AGGT --summary
python analyze.py AGCT AGGT --similarity
python analyze.py AGCT AGGT --edit-distance
python analyze.py AGCT AGGT --edit-distance --load-csv
python analyze.py AGCT AGGT --alignment local
python analyze.py AGCT AGGT --alignment global

python translate.py AUGACGGAGCUUCGGAGCUAG
python translate.py --input-file rna.txt
```

Output examples:
```
python analyze.py ACCC ACCT -e

[[0 1 2 3 4]
 [1 0 1 2 3]
 [2 1 0 1 2]
 [3 2 1 0 1]
 [4 3 2 1 1]]
[['' 'A' 'C' 'C' 'T']
 ['A' '↖' '←' '←' '←']
 ['C' '↑' '↖' '↖' '←']
 ['C' '↑' '↖' '↖' '←']
 ['C' '↑' '↖' '↖' '↖']]
[Edit distance] Cost=1
```
```
python translate.py --input-file rna.txt

MNACFSNLCYESKSIGG
MSDTLSQRLRASLGAIRIAFNLGRSAELD
```
# Customization
Default scoring values: 
```python
# SequenceAnalyzer.py
self.scoring_sys = ScoringSystem(match=1, mismatch=-1, gap=-1)
self.edit_cost_sys = ScoringSystem(match=0, mismatch=1, gap=1)
```
You can set up your own similarity and edit cost matrices by adding `--load-csv` flag:

**scores.csv**

![image](https://user-images.githubusercontent.com/12485656/50089228-7d318480-0205-11e9-9f51-6c396363719d.png)

**edit_cost.csv**

![image](https://user-images.githubusercontent.com/12485656/50089142-43f91480-0205-11e9-8d93-bc05449c039d.png)


