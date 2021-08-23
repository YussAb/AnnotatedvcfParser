# AnnotatedvcfParser.py
Script to parse annotated vcf files (VEP,snpEff,ANNOVAR) to a csv table 

## Documentation for usage from bash terminal 

### Options:

#### Required:
* -v, --vcf    [Path of the annotated VCF file to parse]
* -o, --output [Path of the output CSV table]  

#### Optional:
Split snpEff/VEP transcripts by rows and internal pipe separator by columns
* -t, --splitAnnotator ['s' for snpeff, 'v' for vep, 's,v' for both snpeff and vep together]

Example:
```bash
python3 AnnotatedvcfParser.py -v path/to/AnnotatedVcf.vcf -o path/to/Example.csv  
```

```bash
python3 AnnotatedvcfParser.py -v path/to/VepVcf.vcf -o path/to/Example.csv --splitAnnotator v 
```

```bash
python3 AnnotatedvcfParser.py -v path/to/snpEffVcf.vcf -o path/to/Example.csv --splitAnnotator s 
```

## Documentation for interactive development environment:

You can find all the information inside the jupyter notebook inside this repository. (see:AnnotatedvcfParser.ipynb)

License 
[MIT](https://choosealicense.com/licenses/mit/)

