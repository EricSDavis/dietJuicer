# dietJuicer

## Commands for splitting files

This command splits a 1M read file (4,000,000 lines) into 2 500K read files (using `-l 2000000`):
```
zcat sub_CI_THP1_A_1.1.1_R1.fastq.gz | split -a 3 -l 2000000 -d --additional-suffix=_R1.fastq
zcat sub_CI_THP1_A_1.1.1_R2.fastq.gz | split -a 3 -l 2000000 -d --additional-suffix=_R2.fastq
```

I split the above file into 2 500K read files using above command.

Files were then gzipped with this command:
`for i in $(ls x*); do sbatch -p general --wrap="gzip $i"; done`

Here is an example for how to run external bash scripts:
```
run:
		shell("./scripts/test.sh {input.in1} {output.out1} {log.err1} {input.in2} {output.out2} {log.err2}")
```

```
# slow
num1=$(paste <(gunzip -c x000_R1.fastq.gz) <(gunzip -c x000_R2.fastq.gz) | grep -cE "GATCGATC")
num3=$(paste <(gzcat x000_R1.fastq.gz) <(gzcat x000_R2.fastq.gz) | grep -cE "GATCGATC")

# fast, but inaccurate
num2=$(( $(gzcat x000_R1.fastq.gz | grep -cE "GATCGATC") +  $(gzcat x000_R2.fastq.gz | grep -cE "GATCGATC") ))


paste <(gunzip -c x000_R1.fastq.gz) <(gunzip -c x000_R2.fastq.gz) | grep -cE "GATCGATC"

gunzip -c x000_R1.fastq.gz | grep -n "GATCGATC" | cut -f1 -d : > out1
gunzip -c x000_R2.fastq.gz | grep -n "GATCGATC" | cut -f1 -d : >> out1
num1=$(sort out1 | uniq | wc -l | awk '{print $1}')


gzcat x000_R1.fastq.gz | grep -n "GATCGATC" | cut -f1 -d : > out1
gzcat x000_R2.fastq.gz | grep -n "GATCGATC" | cut -f1 -d : >> out1
num1=$(sort out1 | uniq | wc -l | awk '{print $1}')

num2=$(gunzip -c x000_R1.fastq.gz | wc -l | awk '{print $1}')

```
Final Result:
```
		"""
		gzcat {input.R1} | grep -n {params.ligation} | cut -f1 -d : > {output.temp};
		gzcat {input.R2} | grep -n {params.ligation} | cut -f1 -d : >> {output.temp};
		num1=$(sort {output.temp} | uniq | wc -l | awk '{{print $1}}');
		num2=$(gunzip -c {input.R1} | wc -l | awk '{{print $1}}');

		echo -ne "$num1 " > {output.res};
		echo "$num2" > {output.linecount};
		"""
```

```
awk '{printf("%s %s %s %d %s %s %s %d", $1, $2, $3, 0, $4, $5, $6, 1); for (i=7; i<=NF; i++) {printf(" %s",$i);}printf("\n");}' output/CI_THP1_A_1.1.1_chimera_split1_norm.txt > output/CI_THP1_A_1.1.1_fragment_split1.frag.txt
```
