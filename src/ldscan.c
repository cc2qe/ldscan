#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

double get_chi(double observed,
	       double expected,
	       double min_exp)
{
  double numer;
  if (expected > min_exp && expected != 0) {
    numer = (observed - expected);
    return (numer*numer)/expected;
  }
  else return 0;
}

void get_expected(double *rates_1,
                  double *rates_2,
                  int num_inf,
                  double *d_expected)
{
  int i,j;
  for (i= 0; i < 3; ++i)
    for (j= 0; j < 3; ++j) {
      d_expected[3*i + j] = rates_1[i] * 
	rates_2[j] *
	num_inf;
    }
}

void get_observed(int *locus_1,
		  int *locus_2,
		  int num_samples,
		  double *d_observed,
		  int *multi_informative)
{
  int i;
  for (i = 0; i < 9; ++i) 
    d_observed[i] = 0;

  int genotype;
  for (i = 0; i < num_samples; ++i) {
    // only assess samples that are informative at all k loci
    if (locus_1[i] >= 0 && locus_2[i] >= 0) {
      // multi_informative is the number of samples that are informative
      // at ALL k loci
      *multi_informative += 1;

      // bitwise (tripwise?) representation of genotype
      genotype = 3 * locus_1[i] +
	1 * locus_2[i];
      d_observed[genotype] += 1;
    }
  }
}

void get_rates(int *loci, 
               int num_samples,
               double *rates)
{
  int num_hom_ref = 0,
    num_het = 0,
    num_hom_alt = 0;
  
  int i,loci_i;
  for (i = 0; i < num_samples; ++i) {
    loci_i = loci[i];
    if (loci_i == 0) 
      ++num_hom_ref;
    else if (loci_i == 1)
      ++num_het;
    else if (loci_i == 2)
      ++num_hom_alt;
    // else it's -1 (uninformative)
  }

  int total = num_hom_ref + num_het + num_hom_alt;

  rates[0] = ((double)num_hom_ref) / ((double)total);
  rates[1] = ((double)num_het) / ((double)total);
  rates[2] = ((double)num_hom_alt) / ((double)total);
}

int decToBase(int x,
	      int base)
{
  int i = 0;
  int digit;
  int x_b = 0;

  while(x != 0) {
    digit = x % base;
    x = x / base;

    int tenPower = 1;
    int j;
    for (j = 0; j < i; ++j) {
      tenPower = tenPower * 10;
    }
    x_b += (digit * tenPower);

    ++i;
  }

  return x_b;
}

// number of tests per chunk
// ceil(n*(n-1)/(2c))
unsigned long long int get_chunk_size(unsigned long long int num_loci,
				      int num_chunks)
{
  unsigned long long int num_tests = num_loci * (num_loci - 1) / 2;

  /* fprintf(stderr, "num_tests: %llu\n", num_tests); */

  // this ceilings the number
  unsigned long long int chunk_size;
  chunk_size = (num_tests + num_chunks - 1) / num_chunks;
  return chunk_size;
}

void get_chunk_bounds(unsigned long long int num_loci,
		      int chunk_number,
		      unsigned long long int chunk_size,
		      unsigned long long int *chunk_start,
		      unsigned long long int *chunk_end)
{
  // the preferred number of tests before this chunk
  unsigned long long int pref_tests_start = (chunk_number - 1) * chunk_size;
  unsigned long long int pref_tests_end = chunk_number * chunk_size;

  /* fprintf(stderr, "chunk_size: %llu\n", chunk_size); */
  /* fprintf(stderr, "pref_tests_start: %llu\n", pref_tests_start); */
  /* fprintf(stderr, "pref_tests_end: %llu\n", pref_tests_end); */

  // the i value to start iterating over i,j loci pairs
  unsigned long long int i = 0;
  unsigned long long int tests = 0;
  while (tests < pref_tests_start)
    tests += num_loci - 1 - i++;
  
  // pad each chunk with an extra i to ensure the space is covered
  *chunk_start = i + chunk_number - 1;

  // iterate for the end value
  while (tests < pref_tests_end)
    tests += num_loci - 1 - i++;

  // pad each chunk with an extra i to ensure the space is covered
  *chunk_end = i + chunk_number;

  // control for overflow
  if (*chunk_start > num_loci)
    *chunk_start = num_loci;
  if (*chunk_end > num_loci)
    *chunk_end = num_loci;
  
  return;
  
}

int usage()
{
  fprintf(stderr,
	  "usage: ldscan [options] <file> <samples>\n\n"
	  "author: Ryan Layer and Colby Chiang\n"
	  "version: 0.0.1a\n"
	  "description: finds loci in linkage disequilibrium based on\n"
	  "  genotype frequencies in a set of individuals\n"
	  "\n"
	  "positional arguments:\n"
	  "  file        tab-delimited input file of genotypes\n"
	  "  samples     tab-delimited file of sample names, superpopulations,\n"
	  "                and subpopulations\n"
	  "\n"
	  "optional arguments:\n"
	  "  -h          show this help and exit\n"
	  "  -s INT      number of samples in file\n"
	  "  -l INT      number of loci in file\n"
	  /* "  -k SET_SIZE            number of loci in each set\n" */
	  /* "                           (currently 3 no matter what you put, sucka)\n" */
	  "  -d INT      minimum distance between loci\n"
	  "  -x FLOAT    minimum chi-squared sum to print\n"
	  "  -e FLOAT    cells only contribute to chi-squared value\n"
	  "              if the expected freq is greater than MIN_EXP [5]\n"
	  "  -c INT,INT  process chunk INT out of INT total chunks\n"
	  "  -b          brief but faster output\n"
	  "\n"
	  );
  return 1;
}

int main (int argc, char **argv)
{
  int min_distance = 0;
  uint num_samples;
  unsigned long long int num_loci;
  int set_size;
  double min_chi_sum = 0;
  int yates = 0;
  double min_exp = 5.0;
  int brief = 0;
  
  char raw[100];

  int chunk = 1;
  int num_chunks = 1;
  char *file_name;

  unsigned long long int i;
  int c;
  opterr = 0;

  while ((c = getopt(argc, argv, "hd:s:l:k:x:e:c:b")) != -1) {
    switch (c) {
    case 'h':
      return usage();
    case 'd':
      min_distance = atoi(optarg);
      break;
    case 's':
      num_samples = atoi(optarg);
      break;
    case 'l':
      num_loci = atoi(optarg);
      break;
    case 'k':
      set_size = 2;
      break;
    case 'x':
      min_chi_sum = atoi(optarg);
      break;
    case 'y':
      yates = 1;
      break;
    case 'e':
      min_exp = atof(optarg);
      break;
    case 'b':
      brief = 1;
      break;
    case 'c':
      sscanf(optarg, "%d,%d", &chunk, &num_chunks);
      break;
    case '?':
      if (optopt == 'c')
	fprintf(stderr, "Option -%c requires an argument\n", optopt);
      else if (isprint(optopt))
	fprintf(stderr, "Unknown option '-%c'\n", optopt);
      else
	fprintf(stderr, "Unknown option character '\\x%x'\n", optopt);
      return 1;
    default:
      abort();
    }
  }

  // parse the positional arguments
  file_name = argv[argc-1];
  
  if (argc < 2) {
    return usage();
  }

  // print chunk info
  fprintf(stderr, "Processing chunk %d of %d\n", chunk, num_chunks);

  // the maximum chars per line of file to read in. (should be more than double the number of samples)
  int max_line = 50000;
  char *sep = "\t";
  
  FILE *f = fopen(file_name, "rt");

  // array of arrays containing genotype info for each sample
  // at each locus
  char **chrArr = (char **) malloc(num_loci * sizeof(char*));
  int *posArr = (int *) malloc(num_loci * sizeof(int));
  /* char **geneArr = (char **) malloc(num_loci * sizeof(char*)); */
  char ** rsIdArr = (char **) malloc(num_loci * sizeof(char*));
  int *num_informative = (int *) malloc(num_loci * sizeof(int));
  int **M = (int **) malloc(num_loci * sizeof(int*));
  
  // read each line of the variants file
  // 1       752565  752566  rs3094315       1       2       1       1       2       2
  char line[max_line];
  unsigned long long int j = 0;
  while (fgets(line, max_line, f) != NULL) {
    char *chr = strtok(line, sep);
    int pos = atoi(strtok(NULL, sep)) + 1;
    int end = atoi(strtok(NULL, sep));
    char *rsId = strtok(NULL, sep);
    /* char *gene = strtok(NULL, sep); */
    int inf = 0;
    
    chrArr[j] = strdup(chr);
    posArr[j] = pos;
    /* geneArr[j] = strdup(gene); */
    rsIdArr[j] = strdup(rsId);

    int *locus_gts = (int *) malloc(num_samples * sizeof(int));
    
    i = 0;
    char *tok = strtok(NULL,sep);
    while ((tok != NULL) && (i < num_samples)) {
      locus_gts[i] = atoi(tok);
      tok = strtok(NULL,sep);

      if (locus_gts[i] != -1) {
	++inf;
      }
      ++i;
    }
    
    // store the genotypes into the locus x gt matrix
    M[j] = locus_gts;
    // store the number of informative samples at locus
    num_informative[j] = inf;
    ++j;
  }
  fclose(f);

  // generate array of genotypes (eg 000, 012, 202) so we don't have to
  // call the dec to base function for every locus
  int m_gts[9];
  if (! brief) {
    for (j = 0; j < 9; ++j) {
      m_gts[j] = decToBase(j, 3);
    }
  }

  int k,l;
  double rates_1[3], rates_2[3];
  double expected[9], observed[9], chi[9];
  double chi_sum;
  
  double **rates = malloc(num_loci * sizeof(double));

  /* double *rates[num_loci]; */
  for (i = 0; i < num_loci; ++i) {
    double *rate = (double *) malloc(3 * sizeof(double));
    get_rates(M[i], num_samples, rate);
    rates[i] = rate;
  }

  // get the chunk size, start, and end
  unsigned long long int chunk_size = get_chunk_size(num_loci, num_chunks);
  unsigned long long int chunk_start, chunk_end;

  // chunk_start and chunk_end are updated by this function
  get_chunk_bounds(num_loci,
		   chunk,
		   chunk_size,
		   &chunk_start,
		   &chunk_end);
  
  // analyze the chunk
  fprintf(stderr, "Chunk start locus: %llu\n", chunk_start);
  fprintf(stderr, "Chunk end locus: %llu\n", chunk_end);
  for (i = chunk_start; i < chunk_end; ++i) {
    for (j = i + 1; j < num_loci; ++j) {
      // only calc chi-square if loci are each separated by
      // minimum distance
      if (strcmp(chrArr[i],chrArr[j]) == 0 && abs(posArr[i] - posArr[j]) < min_distance) {
	continue;
      }	
      
      // number of samples at are informative at all loci in k
      int num_multi_informative = 0;
      
      get_observed(M[i],
		   M[j],
		   num_samples,
		   observed,
		   &num_multi_informative);
      
      get_expected(rates[i],
		   rates[j],
		   num_multi_informative,
		   expected);
      
      // calculate chi values for each cell and the chi_sum value for the trio
      chi_sum = 0;
      for (l = 0; l < 9; ++l) {
	chi[l] = get_chi(observed[l],expected[l], min_exp);
	chi_sum += chi[l];
      }
      
      if (chi_sum >= min_chi_sum) {
	if (brief) {
	  printf("%llu\t%llu\t%f\n",
		 i,j,
		 chi_sum);
	}
	
	else {
	  for (l = 0; l < 9; ++l) {
	    printf("%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%s\t%d\t%s\t%.3f\t%.3f\t%.3f\t%02d\t%.0f|%.1f|%.1f\t%f\t%f\n",
		   chrArr[i],posArr[i],rsIdArr[i],rates[i][0],rates[i][1],rates[i][2],
		   chrArr[j],posArr[j],rsIdArr[j],rates[j][0],rates[j][1],rates[j][2],
		   m_gts[l],
		   observed[l],expected[l],observed[l]-expected[l],
		   chi[l],chi_sum);
	  }
	}
      }
    }
  }
  
  for (j = 0; j < num_loci; ++j) {
    free(chrArr[j]);
    /* free(geneArr[j]); */
    free(rsIdArr[j]);
    free(rates[j]);
    free(M[j]);
  }
  free(posArr);
  free(rates);

  return 0;
}
