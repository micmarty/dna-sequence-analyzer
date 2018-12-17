# dna-sequence-analyzer

## Features
- Local alignment (Smith-Waterman algorithm)
- Global alignment (in-progress)
- Edit distance and similarity (Smith-Waterman algorithm)

## Available commands
```
Usage: analyze.py [OPTIONS] SEQUENCE_A SEQUENCE_B

Options:
  -S, --summary
  -s, --similarity
  -e, --edit-distance
  -a, --alignment [global|local]
  --help                          Show this message and exit.
```
## Requirements
- Python 3.7 (type annotations)
- numpy
- click (CLI)

`pip install -r requirements.txt`



## Examples
```
python analyze.py AGCT AGGT --summary
python analyze.py AGCT AGGT --similarity
python analyze.py AGCT AGGT --edit-distance
python analyze.py AGCT AGGT --alignment local
python analyze.py AGCT AGGT --alignment global
```
