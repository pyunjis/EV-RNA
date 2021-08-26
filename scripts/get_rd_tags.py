import pandas as pd

files = snakemake.input
out = snakemake.output[0]

def split_exon_fraction(files, out):
    sample_list = []
    for file in files:
        name = file.split('/')[3]
        results = {}
        equal = False
        for line in open(file):
            if '==' in line:
                equal = True

            if not equal:
                line = line.split()
                key ='_'.join(line[:2])
                results[key] = float(line[-1])

            if equal:
                if not '==' in line and not'Group' in line:
                    line = line.split()
                    key = line[0]
                    results[key] = float(line[2])
                else:
                    continue
        Total_Tags = results['CDS_Exons'] + results["5'UTR_Exons"] + results["3'UTR_Exons"] + results['Introns'] + results['TSS_up_1kb'] + results['TSS_up_5kb'] + results['TSS_up_10kb'] + results['TES_down_1kb'] + results['TES_down_5kb'] + results['TES_down_10kb']
        Exon_Tags = results['CDS_Exons'] + results["5'UTR_Exons"] + results["3'UTR_Exons"]
        Intron_Tags = results['Introns'] 
        Intergenic_Tags = results['TSS_up_1kb'] + results['TSS_up_5kb'] + results['TSS_up_10kb'] + results['TES_down_1kb'] + results['TES_down_5kb'] + results['TES_down_10kb']
        sample_list.append([name, Exon_Tags, Intron_Tags, Intergenic_Tags, Total_Tags])

    frame = pd.DataFrame(sample_list)
    frame.columns = ['Sample','Exon_Tags','Intron_Tags','Intergenic_Tags','Total_Tags']
    frame.to_csv(out,sep='\t',index=False)

split_exon_fraction(files, out)
