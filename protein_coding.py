with open('gencode.v27.onlygrep_protein_coding.gtf') as f:
	with open ('gencode.v27.protein_coding2.gtf','w') as fw:
		for line in f.readlines():
			if line.split('\t')[2]=='gene':
				fw.write(line)
			else:
				if (line.split('\t')[8].split('"')[9])=="protein_coding":
					fw.write(line)

