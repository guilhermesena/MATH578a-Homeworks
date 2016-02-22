#include <stdio.h>
#include <string.h>
#include <algorithm>
using namespace std;

const int ALPHABET_SIZE = 256;
const int MAX_TEXT_SIZE = 3e8;
const int MAX_PATTERN_SIZE = 1e3;
static char text[MAX_TEXT_SIZE];
static char pattern[MAX_PATTERN_SIZE];

static void compute_R(char *pattern, int pattern_length, int *R) {
	for (int i = 0; i < ALPHABET_SIZE; i++) {
		R[i] = pattern_length;
	}
	for (int i = 0; i < pattern_length; i++) {
		R[pattern[i]] = i;
	}
}
static int match(const char *s, int q, int n) {
	int pattern_length = strlen(s);
	for (int i = n; max(q, i) < pattern_length && (s[i] == s[q]); i++, q++)
		;
	return q;
}

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

static void compute_N(int pattern_length, int *N, int *Z) {
	for (int i = 0; i < pattern_length; i++) {
		N[i] = Z[pattern_length - i - 1];
	}
}

static void compute_Lp(int pattern_length, int *Lp, int *N) {
	memset(Lp, -1, sizeof(Lp));
	for (int i = 0; i < pattern_length; i++) {
		int j = pattern_length - N[i] - 1;
		Lp[j] = i;
	}
}

static void compute_lp(int pattern_length, int *lp, int *N) {
	memset(lp, -1, sizeof(lp));
	for (int i = 0; i < pattern_length; i++) {
		int j = N[i];
		lp[j] = i;
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
	compute_R(pattern, pattern_length, R);
	compute_Z(pattern_reverse, pattern_length, Z);
	compute_N(pattern_length, N, Z);
	compute_Lp(pattern_length, Lp, N);
	compute_lp(pattern_length, lp, N);

	int k = pattern_length - 1;
	int m = strlen(text);

	while (k < m) {
		int i, h;
		for (i = pattern_length - 1, h = k; i >= 0 && (pattern[i] == text[h]); i--, h--);
		if (i < 0) {
			printf("Match at [%d, %d]\n", k - pattern_length + 1, k);
			ans++;
			k += pattern_length - lp[1] - 1;
		} else {
			k += max(max(R[text[i]] - i, Lp[i] - i), 1);
		}
	}
	free(R);
	free(Z);
	free(N);
	free(Lp);
	free(lp);

	return ans;
}

int main(int argc, char **argv) {
	scanf(" %s %s",pattern,text);
	printf("MATCHES: %d\n", boyer_moore(pattern,text));
	return 0;
}
