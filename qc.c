#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>
#include <string.h>
#include <pthread.h>
#include <inttypes.h>

#define MAXLINELEN 152

typedef struct {
	char * reads[2];
}fqnum_opt;

typedef struct {
	gzFile inf;
        int64_t total_reads_num;
        int64_t total_base;
        int64_t total_q20;
        int64_t total_q30;
}func_paras;

void use_getopt_long(int argc,char **argv,fqnum_opt * fqs);
void *cal(void * args);
const char * find_file_name(const char *name);

int main(int argc, char * argv[]) {
	fqnum_opt * fqs;
	fqs=(fqnum_opt *)malloc(sizeof(fqnum_opt));
	use_getopt_long(argc,argv,fqs);
	const char * r1_filename=NULL;
	const char * r2_filename=NULL;
	r1_filename=find_file_name(fqs->reads[0]);
	r2_filename=find_file_name(fqs->reads[1]);

	const char * file_read1;
	const char * file_read2;
	time_t start,finish;
   	double totaltime;
	gzFile fp;
	int i;
	func_paras para1,para2;
	pthread_t t1,t2;
	int t1_stat,t2_stat;

	if (((para1.inf)=gzopen(fqs->reads[0],"r"))==NULL){
		fprintf(stderr,"Can not open %s\n",fqs->reads[0]);
		exit(EXIT_FAILURE);
	}
	if (((para2.inf)=gzopen(fqs->reads[1],"r"))==NULL){
		fprintf(stderr,"Can not open %s\n",fqs->reads[1]);
		exit(EXIT_FAILURE);
	}
        para1.total_reads_num=0;
        para1.total_base=0;
        para1.total_q20=0;
        para1.total_q30=0;
        para2.total_reads_num=0;
        para2.total_base=0;
        para2.total_q20=0;
        para2.total_q30=0;
	int64_t sum_reads_num;
	int64_t sum_base;
	int64_t sum_q20;
	int64_t sum_q30;
	int errnum;

	time(&start);
	pthread_create(&t1,NULL,cal,&para1);
	pthread_create(&t2,NULL,cal,&para2);
	t1_stat=pthread_join(t1,NULL);	
	t2_stat=pthread_join(t2,NULL);	
	if(t1_stat!=0){
		fprintf(stderr,"Error: can not join thread1\n");
		exit(EXIT_FAILURE);
	}
	if(t2_stat!=0){
		fprintf(stderr,"Error: can not join thread2\n");
		exit(EXIT_FAILURE);
	}
        if((errnum=gzclose(para1.inf))!=Z_OK){
               fprintf(stderr,"Error when close file %s :%d",fqs->reads[0],errnum);
               exit(EXIT_FAILURE);
        }
        if((errnum=gzclose(para2.inf))!=Z_OK){
               fprintf(stderr,"Error when close file %s :%d",fqs->reads[1],errnum);
               exit(EXIT_FAILURE);
        }


	time(&finish);
	totaltime=(double)finish-start;
	fprintf(stdout,"#FastqName\tReadsNum\tBasesNum\t>=Q20%%\t>=Q30%%\n");
	fprintf(stdout,"%s\t%" PRId64 "\t%" PRId64 "\t%.2f\t%.2f\n",r1_filename,para1.total_reads_num,para1.total_base,(para1.total_q20)*100.0/(para1.total_base),(para1.total_q30)*100.0/(para1.total_base));
	fprintf(stdout,"%s\t%" PRId64 "\t%" PRId64 "\t%.2f\t%.2f\n",r2_filename,para2.total_reads_num,para2.total_base,(para2.total_q20)*100.0/(para2.total_base),(para2.total_q30)*100.0/(para2.total_base));
	sum_reads_num=para1.total_reads_num+para2.total_reads_num;
	sum_base=para1.total_base+para2.total_base;
	sum_q20=para1.total_q20+para2.total_q20;
	sum_q30=para1.total_q30+para2.total_q30;
	fprintf(stdout,"%s\t%" PRId64 "\t%" PRId64 "\t%.2f\t%.2f\n","Total",sum_reads_num,sum_base,sum_q20*100.0/sum_base,sum_q30*100.0/sum_base);
	fprintf(stdout,"##Running Time: %.1fmin\n",totaltime/60);
	return 0;
}

void *cal(void * args) {
	func_paras *p_para;
	p_para=(func_paras *)args;
        char line[MAXLINELEN];
        int tem_base_qual=0;
	int i;
        int64_t tem_total_q20=0;
        int64_t tem_total_q30=0;
	int64_t tem_total_base=0;
	int64_t tem_total_line=0;
	gzFile tem_inf=p_para->inf;
        while (gzgets(tem_inf,line,MAXLINELEN)) {
                (tem_total_line)++;
                if ((tem_total_line)%4!=0)
                        continue;
                else {

                      //  puts(line);
                        i=0;
                        while (line[i]!='\n') {
                                tem_base_qual=line[i]-33;
                                if (tem_base_qual>=30){
                                        tem_total_q30++;
                                        tem_total_q20++;
                                }
                                else if (tem_base_qual>=20)
                                        tem_total_q20++;
                                (tem_total_base)++;
                                i++;
                        }
                }
        }
        (p_para->total_q20)=tem_total_q20;
        (p_para->total_q30)=tem_total_q30;
	(p_para->total_reads_num)=tem_total_line/4;
	(p_para->total_base)=tem_total_base;
}

void use_getopt_long(int argc,char **argv,fqnum_opt * fqs){

	const char *optstring = "1:2:hv";
	int c;
	struct option opts[] = {
		{"read1", 1, NULL, '1'},
		{"read2", 1, NULL, '2'},
		{"version",0,NULL,'v'},
		{"help", 0, NULL, 'h'}
	};
	while((c = getopt_long(argc,argv,optstring,opts,NULL)) != -1){
		switch(c){
			case '1':
				fqs->reads[0]=optarg;
				break;
			case '2':
				fqs->reads[1]=optarg;
				break;
			case 'v':
				fprintf(stdout,"version-20190529 \n");
				exit(EXIT_FAILURE);
			case 'h':
				fprintf(stdout,"\n");
				fprintf(stdout,"Usage: FastqNums -1 <read1> -2 <read2>\n");
				fprintf(stdout,"        -1	--read1	STR	read1 path\n");
				fprintf(stdout,"        -2	--read2	STR	read2 path\n");
				fprintf(stdout,"        -v	--version	software version\n");
				fprintf(stdout,"        -h	--help		help message\n");
				exit(EXIT_SUCCESS);
			case '?':
				fprintf(stderr,"Unknown option, you may need \'-h\' or \'--help\'\n");
				exit(EXIT_FAILURE);
			default:
				fprintf(stderr,"Unknown option, you may need \'-h\' or \'--help\'\n");
				exit(EXIT_FAILURE);	
		}
	}
	if(NULL==fqs->reads[0] && NULL==fqs->reads[1]){
		fprintf(stdout,"\n");
                fprintf(stdout,"Usage: FastqNums -1 <read1> -2 <read2>\n");
                fprintf(stdout,"        -1      --read1 STR     read1 path\n");
                fprintf(stdout,"        -2      --read2 STR     read2 path\n");
                fprintf(stdout,"        -v      --version       software version\n");
                fprintf(stdout,"        -h      --help          help message\n");
                exit(EXIT_SUCCESS);
	}
	if(NULL==fqs->reads[0]){
		fprintf(stderr,"\n");
		fprintf(stderr,"Unknown option, you may need \'-h\' or \'--help\'\n");
		exit(EXIT_FAILURE);	
	}
	if(NULL==fqs->reads[1]){
		fprintf(stderr,"\n");
		fprintf(stderr,"Unknown option, you may need \'-h\' or \'--help\'\n");
		exit(EXIT_FAILURE);	
	}
}
const char * find_file_name(const char *name){
	const char *name_start = NULL;
	int sep = '/';
	if (NULL == name) {
		fprintf(stderr,"the path name is NULL\n");
		return NULL;
	}	
	name_start = strrchr(name, sep);
	return (NULL == name_start)?name:(name_start + 1);
}
