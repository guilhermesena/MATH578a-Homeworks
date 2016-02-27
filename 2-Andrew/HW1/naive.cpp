#include <stdio.h>
#include <string.h>
const int MAX_TEXT_SIZE = 3e8;
const int MAX_PATTERN_SIZE = 1e3;
static char P[MAX_PATTERN_SIZE];
static char T[MAX_TEXT_SIZE];
long long int num_comparisons = 0;
int naive(char *pattern, char *text) {
	int lp = strlen(pattern), lt = strlen(text);
	int ans = 0;
	for (int i = 0; i < lt - lp + 1; i++) {
		for (int j = 0; j < lp; j++) {
			num_comparisons++;
			if (pattern[j] != text[i + j])
				break;
			if (j == lp - 1) {
				ans++;
				printf("Match at [%d, %d]\n", i, i + lp - 1);
			}
		}
	}
	return ans;
}

int main(int argc, char **argv) {
	if (argc != 2) {
		printf("Usage: ./naive PATTERN <chromosome_input");
		return 1;
	}
	strcpy(P, argv[1]);
	scanf(" %s", T);
	printf("MATCHES: %d\n", naive(P, T));
	printf("Num comparisons: %lld\n", num_comparisons);
	return 0;
}
