from ..version import __version__

def main():
    print("Welcome to CamoTSS v%s! Command lines available:\n " %(__version__))
    print("CamoTSS\n    Counting reads for each transcript from "
          "per cell bam files")

if __name__ == "__main__":
    main()