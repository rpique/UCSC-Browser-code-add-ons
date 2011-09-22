#ifndef PWM_H
#define PWM_H

#define maxPSSM 200
#define maxStr 500

struct pssm
{
  void *next;
  char *name;
  int w;
  double **matrix; /* */ 
};

void readInJasparPwm(struct pssm *pm, char *fileName);

void readInTransfacPwm(struct pssm *pm, char *fileName);
void allocateMemoryToMatrix (struct pssm *pwm);
void convertPSSMToLogs (struct pssm *pwm);
void reversePWM(struct pssm *op,const struct pssm *p);

//void initialise_pssm(struct pssm *pm,char *fileName);
void initialise_pssm(struct pssm *pwm,char *fileMotif,double addPseudoCounts, boolean useJaspar);


void printMatrix (struct pssm *pm);  
//void read_seq_fasta(struct parameters *p, struct sequence *s);
//void print_seq_fasta(FILE *f, struct sequence s);
int ATGCbase(char b, boolean maskRep);
double compare_subseq_to_pssm (char *seq, struct pssm *p, boolean useSnpRobust);

//int get_pwm_scores(struct seq *seq, struct pssm *p);



#endif /* PWM_H */

