#define _POSIX_C_SOURCE 199309L

#include "GraphBLAS.h"
#undef I

#include <inttypes.h>
#include <malloc.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>

#define TRY_GrB(method, rc, msg, label) { \
	GrB_Info info = method; \
	if (!(info == GrB_SUCCESS || info == GrB_NO_VALUE)) { \
		printf ("Error: %s\n", msg); \
		printf ("    GraphBLAS error: %d\n", info); \
		printf ("    Failure: line %d file %s\n", __LINE__, __FILE__); \
		rc = EXIT_FAILURE; \
		goto label; \
	} \
} \

#define TRY(method, rc, msg, label) { \
	int info = method; \
	if (info != EXIT_SUCCESS) { \
		printf ("Error: %s\n", msg); \
		printf ("    Failure: line %d file %s\n", __LINE__, __FILE__); \
		rc = EXIT_FAILURE; \
		goto label; \
	} \
} \

#define EL_TYPE(SIZE) uint ## SIZE ## _t
#define DEFINE_KRON_OP(SIZE) \
void kron_op_ ## SIZE(void *z, const void *x, const void *y) { \
	const EL_TYPE(SIZE) n_mask = (EL_TYPE(SIZE))((EL_TYPE(SIZE))0xe0 << (SIZE - 8)); \
	const EL_TYPE(SIZE) t_mask = (EL_TYPE(SIZE))(~n_mask); \
	const EL_TYPE(SIZE) a = * (EL_TYPE(SIZE) *) x; \
	const EL_TYPE(SIZE) b = * (EL_TYPE(SIZE) *) y; \
	* (bool *) z = (n_mask & a & b) \
		|| ((t_mask & a) && ((t_mask & a) == (t_mask & b))); \
} \

#define KRON_OP(SIZE) kron_op_ ## SIZE

DEFINE_KRON_OP(8)
DEFINE_KRON_OP(16)
DEFINE_KRON_OP(32)
DEFINE_KRON_OP(64)

static void print_usage(const char *name) {
	fprintf(stderr, "Usage %s <RSA> <G> <M_meta> <Out>\n", name);
	fprintf(stderr, "\n"
	"- [in] RSA -- path to csv file (src,dst,val) containing matrix A\n"
	"- [in] G -- path to csv file (src,dst,val) containing matrix B\n"
	"- [in] M_meta -- path to text file containing data about matrices\n"
	"    - n_bits -- number of nonterminal bits\n"
	"    - t_bits -- number of terminal bits: n_bits + t_bits == 8|16|32|64\n"
	"    - G_n -- size of square matrix G\n"
	"    - RSA_n -- size of square matrix RSA\n"
	"- [in] E_metadata -- path to file with experiment metadata"
	"    - Warmup reps -- number of warmup runs\n"
	"    - Measure reps -- number of measurement runs\n"
	"- [out] Out -- path to file for saving benchmark results\n"
	 );
}

struct kron_input {
	GrB_Type input_elemets_type;
	void *operation;
};

static int get_kron_input(struct kron_input *k, int bits) {
	switch (bits) {
		case 8:
			k->input_elemets_type = GrB_UINT8;
			k->operation = KRON_OP(8);
			break;
		case 16:
			k->input_elemets_type = GrB_UINT16;
			k->operation = KRON_OP(16);
			break;
		case 32:
			k->input_elemets_type = GrB_UINT32;
			k->operation = KRON_OP(32);
			break;
		case 64:
			k->input_elemets_type = GrB_UINT64;
			k->operation = KRON_OP(64);
			break;
		default:
			return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

struct matr_metadata {
	int bits;
	size_t rsa_size;
	size_t graph_size;
	size_t res_size;
};

struct exp_metadata {
	size_t w_runs;
	size_t m_runs;
};

static int read_m_metadata(struct matr_metadata *m, const char *filename) {
	FILE *f = fopen(filename, "r");
	if (!f) {
		return EXIT_FAILURE;
	}
	int rc = EXIT_SUCCESS;
	if (fscanf(f, "%d", &m->bits) != 1) {
		rc = EXIT_FAILURE;
		goto read_m_metadata_out;
	}
	int t_bits;
	if (fscanf(f, "%d", &t_bits) != 1) {
		rc = EXIT_FAILURE;
		goto read_m_metadata_out;
	}
	m->bits += t_bits;

	if (fscanf(f, "%zu", &m->graph_size) != 1) {
		rc = EXIT_FAILURE;
		goto read_m_metadata_out;
	}
	if (fscanf(f, "%zu", &m->rsa_size) != 1) {
		rc = EXIT_FAILURE;
		goto read_m_metadata_out;
	}
	m->res_size = m->rsa_size * m->graph_size;
read_m_metadata_out:
	fclose(f);
	return rc;
}

static int read_e_metadata(struct exp_metadata *m, const char *filename) {
	FILE *f = fopen(filename, "r");
	if (!f) {
		return EXIT_FAILURE;
	}
	int rc = EXIT_SUCCESS;
	if (fscanf(f, "%zu", &m->w_runs) != 1) {
		rc = EXIT_FAILURE;
		goto read_m_metadata_out;
	}
	if (fscanf(f, "%zu", &m->m_runs) != 1) {
		rc = EXIT_FAILURE;
		goto read_m_metadata_out;
	}
read_m_metadata_out:
	fclose(f);
	return rc;
}

static int read_matrix(GrB_Matrix *M, const char *filename, size_t ndim, GrB_Type type) {
	int rc = 0;
	FILE *f = fopen(filename, "r");
	if (!f) {
		fprintf(stderr, "Cannot open file '%s'\n", filename);
		return 1;
	}
	size_t len = 256;
	size_t ntuples = 0;
	size_t xsize;
	GxB_Type_size(&xsize, type);
	GrB_Index *I = malloc(len * sizeof(I[0]));
	GrB_Index *J = malloc(len * sizeof(J[0]));
	void *X = malloc (len * xsize);
	if (!I || !J || !X) {
		fprintf(stderr, "Cannot allocate memory for csv entries");
		rc = EXIT_FAILURE;
		goto read_matrix_free_vectors;
	}
	GrB_Index *I2, *J2;
	void *X2;
	GrB_Index i, j;
	uint64_t val;
	while (fscanf(f, "%"SCNd64",%"SCNd64",%"SCNu64, &i, &j, &val) != EOF) {
		if (ntuples >= len) {
			I2 = realloc(I, 2 * len * sizeof(I[0]));
			J2 = realloc(J, 2 * len * sizeof(J[0]));
			X2 = realloc(X, 2 * len * xsize);
			if (I2 == NULL || J2 == NULL || X2 == NULL) {
				fprintf(stderr, "Cannot reallocate memory for csv entries");
				rc = EXIT_FAILURE;
				goto read_matrix_free_vectors;
			}
			I = I2; I2 = NULL;
			J = J2; J2 = NULL;
			X = X2; X2 = NULL;
			len = len * 2 ;
		}
		I[ntuples] = i;
		J[ntuples] = j;
		memcpy(X + ntuples * xsize, &val, xsize);
		ntuples++;
	}
	if (GrB_Matrix_new(M, type, ndim, ndim) != GrB_SUCCESS) {
		fprintf(stderr, "Cannot create Matrix object\n");
		rc = EXIT_FAILURE;
		goto read_matrix_free_vectors;
	}
	if (GrB_Matrix_build(*M, I, J, X, ntuples, GxB_IGNORE_DUP) != GrB_SUCCESS) {
		fprintf(stderr, "Cannot build Matrix from vectors\n");
		rc = EXIT_FAILURE;
		goto read_matrix_free_vectors;
	}
read_matrix_free_vectors:
	free(I);
	free(J);
	free(X);
read_matrix_out:
	fclose(f);
	return rc;
}

static int mem_get_hwm(size_t *m) {
	malloc_trim(0);
	FILE *f = fopen("/proc/self/status", "r");
	if (!f) {
		fprintf(stderr, "Cannot read /proc/self/status\n");
		return EXIT_FAILURE;
	}
	char line[256];
	while (fgets(line, sizeof(line), f)) {
		if (sscanf(line, "VmHWM: %zu kB", m) == 1) {
			break;
		}
	}
	fclose(f);
	return EXIT_SUCCESS;
}

static int mem_get_rss(size_t *m) {
	malloc_trim(0);
	FILE *f = fopen("/proc/self/status", "r");
	if (!f) {
		fprintf(stderr, "Cannot read /proc/self/status\n");
		return EXIT_FAILURE;
	}
	char line[256];
	while (fgets(line, sizeof(line), f)) {
		if (sscanf(line, "VmRSS: %zu kB", m) == 1) {
			break;
		}
	}
	fclose(f);
	return EXIT_SUCCESS;
}

static int mem_reset_hwm() {
	malloc_trim(0);
	FILE *f = fopen("/proc/self/clear_refs", "w");
	if (!f || fputs("5\n", f) == EOF) {
		fprintf(stderr, "Cannot write /proc/self/clear_refs\n");
		return EXIT_FAILURE;
	}
	fclose(f);
	return EXIT_SUCCESS;
}

#define MAX_RUNS 100

struct experiment_result {
	double time;
	size_t mem_before;
	size_t mem_peak;
	size_t mem_after;
	size_t nnz;
};

static struct matr_metadata m_metadata;
static struct exp_metadata e_metadata;

static GrB_Matrix Rsa = NULL;
static GrB_Matrix Graph = NULL;
static GrB_BinaryOp kron_binop = NULL;
static size_t csize;

/* селектор -- побайтовая проверка на 0
первый аргемент -- адрес, куда нужно записать bool
второй аргумент -- адрес значения, которое хотим фильтровать*/
void simple_fun(void *z, const void *x)
{
	bool result = false;

	const unsigned char *bytes = (const unsigned char *)x;
	for (size_t i = 0 ; i < csize ; ++i)
	{
		if (*(bytes + i))
		{
			result = true;
			break;
		}
	}

	*(bool *)z = result;
}
// приводим пользовательскую функцию к стандартной GxB_unary_function
static GxB_unary_function selector = (GxB_unary_function)simple_fun;

static struct experiment_result results[MAX_RUNS];

static int run_kron(size_t iterations) {
	int rc = EXIT_SUCCESS;
	GrB_Matrix C = NULL;
	struct timespec t_start, t_end;
	for (size_t it = 0; it < iterations; ++it) {
		TRY_GrB(GrB_Matrix_new (&C, GrB_BOOL, m_metadata.res_size, m_metadata.res_size), rc, "Cannot create matrix for result\n", run_kron_out);
		TRY(mem_get_rss(&results[it].mem_before), rc, "Cannot get RSS before kron\n", run_kron_mat_free);
		TRY(mem_reset_hwm(), rc, "Cannot reset HWM before kron\n", run_kron_mat_free);

		// определяем размер типа для работы селектора
		GrB_Type ctype;
		GxB_BinaryOp_ztype(&ctype, kron_binop); 
		size_t csize_cur;
		GxB_Type_size(&csize_cur, ctype);
		csize = (size_t)csize_cur;

		clock_gettime(CLOCK_MONOTONIC, &t_start);
		TRY_GrB(GrB_Matrix_kronecker_BinaryOp (C, NULL, NULL, kron_binop, selector, Rsa, Graph, NULL), rc, "Cannot calculate kron\n", run_kron_out);
		TRY_GrB(GrB_Matrix_nvals (&results[it].nnz, C), rc, "Cannot get nvals after kron\n", run_kron_out);
		clock_gettime(CLOCK_MONOTONIC, &t_end);
		/*if (it == 0) {
			GrB_Index cnvals = results[it].nnz;
			GrB_Index *I_out = malloc(cnvals * sizeof(GrB_Index));
			GrB_Index *J_out = malloc(cnvals * sizeof(GrB_Index));
			double   *X_out = malloc(cnvals * sizeof(double));
			if (I_out && J_out && X_out) {
				GrB_Matrix_extractTuples_FP64(I_out, J_out, X_out, &cnvals, C);
				FILE *fout = fopen("result_C.csv", "w");
				if (fout) {
					for (GrB_Index k = 0; k < cnvals; k++)
						fprintf(fout, "%" PRIu64 ",%" PRIu64 "\n", I_out[k], J_out[k]);
					fclose(fout);
				}
				free(I_out); free(J_out); free(X_out);
			}
		}*/
		TRY(mem_get_hwm(&results[it].mem_peak), rc, "Cannot get RSS after kron\n", run_kron_mat_free);

		results[it].time = (t_end.tv_sec - t_start.tv_sec) + (t_end.tv_nsec - t_start.tv_nsec) / 1e9;
		TRY(mem_get_rss(&results[it].mem_after), rc, "Cannot get RSS after kron\n", run_kron_mat_free);
		TRY_GrB(GrB_Matrix_free(&C), rc, "Cannot clean result matrix\n", run_kron_out);
	}
run_kron_out:
	return rc;
run_kron_mat_free:
	GrB_Matrix_free(&C);
	return rc;
}

static int write_results(const char *filename) {
	FILE *f = fopen(filename, "w");
	if (!f) {
		fprintf(stderr, "Cannot open file '%s' for writing results\n", filename);
		return EXIT_FAILURE;
	}
	fputs("iteration,nnz,time,rss_before,hwm,rss_after\n", f);
	for (size_t i = 0; i < e_metadata.m_runs; ++i) {
		fprintf(f, "%zu,%zu,%.6f,%zu,%zu,%zu\n", i, results[i].nnz, results[i].time, results[i].mem_before, results[i].mem_peak, results[i].mem_after);
	}
	fclose(f);
	return EXIT_SUCCESS;
}

int main(int argc, char *argv[]) {
	int rc = EXIT_SUCCESS;
	if (argc != 6) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
	TRY(read_m_metadata(&m_metadata, argv[3]), rc, "Cannot get matrices metadata from file\n", out);
	TRY(read_e_metadata(&e_metadata, argv[4]), rc, "Cannot get experiment metadata from file\n", out);
	struct kron_input kron_settings;
	TRY(get_kron_input(&kron_settings, m_metadata.bits), rc, "Unsupported number of bits\n", out);

	TRY_GrB(GrB_init(GrB_NONBLOCKING), rc, "Cannot initialize GraphBLAS\n", out);

	GrB_Global_set_INT32 (GrB_GLOBAL, true, GxB_BURBLE) ;
	TRY_GrB(GrB_BinaryOp_new(&kron_binop, kron_settings.operation, GrB_BOOL, kron_settings.input_elemets_type, kron_settings.input_elemets_type), rc, "Cannot create binary operation\n", out_fin_graphblas) ;

	TRY(read_matrix(&Rsa, argv[1], m_metadata.rsa_size, kron_settings.input_elemets_type), rc, "Cannot read RSA\n", out_free_binop);
	TRY(read_matrix(&Graph, argv[2], m_metadata.graph_size, kron_settings.input_elemets_type), rc, "Cannot read Graph\n", out_free_rsa);
	// Warmup
	TRY(run_kron(e_metadata.w_runs), rc, "Cannot run warmup iterations\n", out_free_graph);
	// Measure
	TRY(run_kron(e_metadata.m_runs), rc, "Cannot run measure iterations\n", out_free_graph);
	// Write results
	TRY(write_results(argv[5]), rc, "Cannot write results to file\n", out_free_graph);
out_free_graph:
	GrB_Matrix_free(&Graph);
out_free_rsa:
	GrB_Matrix_free(&Rsa);
out_free_binop:
	GrB_BinaryOp_free(&kron_binop);
out_fin_graphblas:
	GrB_finalize();
out:
	return rc;
}
