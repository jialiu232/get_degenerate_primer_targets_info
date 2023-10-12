import argparse

parser = argparse.ArgumentParser(description="Summarize the mfpqpcr evaluation result for single primers to degenerate primers.")
parser.add_argument("-i", "--input", type=str, help="input file from mfpqpcr results")
parser.add_argument("-fo", "--output_1", action='store_true', help="output file 1 name")
parser.add_argument("-ft", "--output_2", action='store_true', help="output file 2 name")
args = parser.parse_args()


## Define the first output file argument
#print(args.output_1)
fout_one = ""
if args.output_1:
    fout_one = args.output_1
else:
    fout_one = f'{args.input}.withPrimerCluster.txt'
#print(fout_one)

## Define the second output file argument
fout_two = ""
if args.output_2:
    fout_two = args.output_2
else:
    fout_two = f'{args.input}.ranked_primerCluster_tarCounts.txt'
#print(fout_two)



## a function to extract the degenerate primer cluster name
def extract_primerClus(lin):
    lineList = lin.split()
    prims = f'{lineList[0]},{lineList[2]}'
    #print(prims)

    primClus = [p.split('.')[-2] for p in prims.split(',')]
    #print(primClus)

    finalCluster = [x for x in set(primClus) if primClus.count(x) > 1][0]
    return(finalCluster)


## summarize the degenerate primer pairs for the evaluated single primer pairs

nameL = ""
clus = ""
prim_tar_dic = {}

out1 = open(fout_one, 'w')

with open(args.input, 'r') as fin:
    for line in fin:
        line = line.strip()
        if line[0] == 'F':
            #print(line)
            clus = extract_primerClus(line)
            nameL = line
            #print(clus)
        elif line[0] == '>':
            #print(f'{clus}\t{nameL}\t{line}')
            out1.write(f'{clus}\t{nameL}\t{line}\n')

            if clus not in prim_tar_dic:
                prim_tar_dic[clus] = [line]
            else:
                prim_tar_dic[clus].append(line)


out1.close()
#print(prim_tar_dic)


## Write the sorted best primers and the number of targets into a file
prim_tarCounts = {key:len(val) for (key,val) in prim_tar_dic.items()}
#print(prim_tarCounts)

N = len(prim_tarCounts)
sorted_prim_tarCounts = dict(sorted(prim_tarCounts.items(), key = lambda x: x[1], reverse = True)[:N])
#print(sorted_prim_tarCounts)

out2 = open(fout_two, 'w')

for key in sorted_prim_tarCounts:
    out2.write(f'{key}\t{sorted_prim_tarCounts[key]}\n')

out2.close()
