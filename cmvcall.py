from script.run import Run
from pathlib import Path

def main():
    classrun = Run()
    args = classrun.get_args()
    inputdir = Path(args.input)
    cinput = args.cinput
    correct = args.CORRECT
    a1 = float(args.a1)
    a2 = float(args.a2)
    FC = float(args.FC)
    depth = float(args.depth)
    classrun.run(inputdir,correct,a1,a2,FC,depth,cinput)
if __name__=="__main__":
    main()
