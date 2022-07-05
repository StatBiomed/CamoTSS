from ..version import __version__

def main():
    print("Welcome to scTSS v%s! Command lines available:\n " %(__version__))
    print("scTSS-count\n    Counting reads for each transcript from "
          "per cell bam files")
    print("scTSS-quant\n    Quantify alternative transcription start site ratio" 
           "and detecting differential transcription start site\n")

if __name__ == "__main__":
    main()