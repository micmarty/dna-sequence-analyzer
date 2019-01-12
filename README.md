![](https://img.shields.io/badge/python-3.7-blue.svg?style=popout-square)
![](https://img.shields.io/badge/platform-Linux_|%20Windows%20|_macOS-blue.svg?style=popout-square)
![](https://img.shields.io/badge/license-Apache%202.0-blue.svg?style=popout-square)

# Sequence analyzer (bioinformatics)

## About
This is a set of simple command-line python scripts with [101](https://dictionary.cambridge.org/dictionary/english/101) algorihtms used in bioinformatics.

## Features
- Pariwise local alignment (Smith-Waterman algorithm)
- Pairwise global alignment (Needleman-Wunsch (with or without start gap penalty) and Hirschberg's algorithm)
- Edit distance and similarity (Needleman-Wunsch algorithm)
- RNA to amino acids translation

## Available commands
```
Usage: analyze.py [OPTIONS] SEQUENCE_A SEQUENCE_B

Options
  -S, --summary
  -s, --similarity
  -e, --edit-distance
  -a, --alignment [global|local]
  --gap-start INTEGER             When gap starts, that score will be applied
                                  (besides `gap` penalty).
  --gap INTEGER                   Penalty for a single gap.
  --match INTEGER                 Score for identical letters in both
                                  seqences.
  --mismatch INTEGER              Score for mismatched letters in both
                                  seqences.
  --load-csv                      Load scores.csv and edit_cost.csv
  --verbose                       Verbose mode.
  --help                          Show this message and exit.
```
```
Usage: translate.py [OPTIONS] [SEQUENCE]

Options:
  -i, --input-file FILE  Path to text file containing long nucleotide sequences (1 sequence = 1 line)
  --help                 Show this message and exit.
```

## Usage examples
```
python analyze.py AGCT AGGT --summary
python analyze.py AGCT AGGT --similarity
python analyze.py AGCT AGGT --edit-distance
python analyze.py AGCT AGGT --edit-distance --load-csv
python analyze.py AGCT AGGT --alignment local
python analyze.py AGCT AGGT --alignment global 

# TODO Test these in different configuration
python analyze.py AGCT AGGT --alignment global --gap-start -10 --gap -1
python analyze.py AGCT AGGT --alignment global --match 10 --mismatch -10
python analyze.py AGCT AGGT --alignment global --load_csv --gap-start -3

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

## Requirements
- Python 3.7 (type annotations)
- numpy (storing matrices)
- pandas (loading CSV into DataFrame)
- click (CLI interface)

We recommend using `conda`/`virtualenv`/`pyenv` environment (this step is optional)

`conda create --name sequence-analyzer-env python=3.7 pip`

## Requirements installation

`pip install -r requirements.txt`


## Customization
Default scoring values: 
```python
# SequenceAnalyzer.py
self.scoring_sys = ScoringSystem(match=1, mismatch=-1, gap=-1)
self.edit_cost_sys = ScoringSystem(match=0, mismatch=1, gap=1)
```
You can set up your own similarity and edit cost matrices by adding `--load-csv` flag

(these files are read by default)

**scores.csv**

![image](https://user-images.githubusercontent.com/12485656/50089228-7d318480-0205-11e9-9f51-6c396363719d.png)

**edit_cost.csv**

![image](https://user-images.githubusercontent.com/12485656/50089142-43f91480-0205-11e9-8d93-bc05449c039d.png)

(Note: if any of your sequences contains invalid symbols, default values from `ScoringSystem` will be used instead)

## Credits
- ![](https://avatars2.githubusercontent.com/u/12485656?s=22&v=4) [Michał Martyniak (@micmarty)](http://martyniak.me)
- Artur Śliwa [(@asliwa)](https://github.com/asliwa)

## License

Feel free play around with our code. If you see any bugs, please tell us about them in issues :heart:!

[Apache License 2.0](http://www.apache.org/licenses/LICENSE-2.0.html)
