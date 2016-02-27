#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;

const int ALPHABET_SIZE = 256;
const int MAX_TEXT_SIZE = 3e8;
const int MAX_PATTERN_SIZE = 1e3;
static char text[MAX_TEXT_SIZE];
static char pattern[MAX_PATTERN_SIZE];
long long int num_comparisons = 0;

//Method used to reverse a string and calculate Z
static void strrev(char *p) {
	char *q = p;
	while (q && *q)
		++q;
	for (--q; p < q; ++p, --q)
		*p = *p ^ *q, *q = *p ^ *q, *p = *p ^ *q;
}

//R is the rightmost position of each ASCII character in the pattern. 0 if not present
static void compute_R(char *pattern, int pattern_length, int *R) {
	memset(R, 0, sizeof(R));
	for (int i = 0; i < pattern_length; i++) {
		R[pattern[i]] = i;
	}
}

//Match function for the Z computation
static int match(const char *s, int q, int n) {
	int pattern_length = strlen(s);
	for (int i = n; max(q, i) < pattern_length && (s[i] == s[q]); i++, q++)
		;
	return q;
}

//Z algorithm mostly taken form the cpp provided. Used to calculate N
static void compute_Z(char *pattern, int pattern_length, int *Z) {
	int l = 0, r = 0;
	for (int k = 1; k < pattern_length; k++) {
		if (k >= r) {
			Z[k] = match(pattern, 0, k);
			if (Z[k] > 0) {
				r = k + Z[k];
				l = k;
			}
		} else {
			if (Z[k - l] < r - k) {
				Z[k] = Z[k - l];
			} else {
				int q = match(pattern, r, r - k);
				Z[k] = q - k;
				r = q;
				l = k;
			}
		}
	}
}

//N[j] is the length of the longest suffix of the substring P[1..j] that
//Is also a suffix of P. Used to calculate Lp and lp
static void compute_N(int pattern_length, int *N, int *Z) {
	for (int i = 0; i < pattern_length; i++) {
		N[i] = Z[pattern_length - i - 1];
	}
}

//Lp[i] is the largest index j less than |P| such that N[j] = |P[i..n]|
static void compute_Lp(int pattern_length, int *Lp, int *N) {
	memset(Lp, 0, sizeof(Lp));
	for (int i = 0; i < pattern_length; i++) {
		int j = pattern_length - N[i] - 1;
		Lp[j] = i;
	}
}

//lp(i) is the largest j <= n - i + 1 such that N[j] = j
static void compute_lp(int pattern_length, int *lp, int *N) {
	memset(lp, 0, sizeof(lp));
	int i = 0, j = pattern_length - 1;
	while (j >= 0 && i < pattern_length) {
		if (N[j] == j + 1) {
			while (j <= pattern_length - i + 1 && i < pattern_length) {
				lp[i++] = j + 1;
			}
		}
		j--;
	}
}

static int boyer_moore(char *pattern, char *text) {
	int ans = 0;
	const int pattern_length = strlen(pattern);
	int *R = new int[ALPHABET_SIZE];
	int *Z = new int[pattern_length];
	int *N = new int[pattern_length];
	int *Lp = new int[pattern_length];
	int *lp = new int[pattern_length];
	char pattern_reverse[pattern_length + 1];
	strcpy(pattern_reverse, pattern);
	strrev(pattern_reverse);
	compute_R(pattern, pattern_length, R);
	compute_Z(pattern_reverse, pattern_length, Z);
	compute_N(pattern_length, N, Z);
	compute_Lp(pattern_length, Lp, N);
	compute_lp(pattern_length, lp, N);
	int k = pattern_length - 1;
	int m = strlen(text);

	while (k < m) {
		int i = pattern_length - 1, h = k;
		for (; (i >= 0) && (pattern[i] == text[h]); i--, h--)
			num_comparisons++;
		if (i < 0) { //Match found
			printf("Match at [%d, %d]\n", k - pattern_length + 1, k);
			ans++;
			k += pattern_length - lp[1];
		} else {
			//Skip as much as possible between bad prefix character and good suffix
			k += max(max(1, i - R[text[h]]),
					pattern_length - 1 - ((Lp[i] > 0) ? (Lp[i]) : (lp[i])));
		}
	}
	delete [] R
	delete [] Z
	delete [] N
	delete [] Lp
	delete [] lp
	return ans;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		printf("Usage: ./boyer-moore PATTERN <chromosome_file");
		return 1;
	}
	strcpy(pattern, argv[1]);
	scanf(" %s", text);
	printf("Pattern size = %d, Text size = %d\n", (int) strlen(pattern),
			(int) strlen(text));
	printf("MATCHES: %d\n", boyer_moore(pattern, text));
	printf("Num comparisons: %lld\n", num_comparisons);
	return 0;
}
