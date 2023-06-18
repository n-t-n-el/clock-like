# clock-like
Age-dependence of mutation in germline and soma  


### Datasets

### ```data/germline.csv.gz```

Total of 714267 SNVs. Partitioned by original datasets

```
An2018             238163
Francioli2014       10773
Goes2019             6363
Goldmann2016        35793
Halldorsson2019    180154 (includes Jonsson2017)
Li2017              35471
Lin2021              3287
Michaelson2012        581
Rahbari2015           747
Richter2020         54011
Sasani2019          28844
Yuen2017           120080
```

```
An2018              https://www.science.org/doi/10.1126/science.aat6576
Francioli2014        
Goes2019              
Goldmann2016        https://www.nature.com/articles/s41588-018-0071-6 
Halldorsson2019     https://www.science.org/doi/10.1126/science.aau1043
Jonsson2017				  
Li2017              
Lin2021              
Michaelson2012        
Rahbari2015           
Richter2020         
Sasani2019          https://elifesciences.org/articles/46922
Yuen2017           
```

An2018, Goldmann2016, Halldorsson2019 compiled by Ipsita Agarwal (An2018 and Goldmann2016 phasing: private communication)
All unspecified downloaded via Gene4Denovo website

```
Zhao2019            https://doi.org/10.1093/nar/gkz923
```

Content: column names

```
Chr                   chromosome
Pos                   position (hg19)
Ref       
Alt				
Context               3-mer reference context
Type                  one of 96 mutation types
Phase                 mother/father
PhasingProba          probability of phasing in given proband 
Dataset							  
Proband               offspring ID
Father                father ID
Mother                mother ID
FatherAgeAtBirth	    
MotherAgeAtBirth
Phenotype             case/control/unascertained 
```

All unspecified should be selbstverstaendig 
All phased mutations are unascertained. There is a total of 

```
123375     paternal mutations
33434 	   maternal mutations
```

Probability of phasing in given proband was estimated as a fraction of phased de novo mutations.

