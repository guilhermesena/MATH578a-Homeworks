#include <stdio.h>
#include <string.h>

int naive(char *pattern, char *text) {
	int lp = strlen(pattern), lt = strlen(text);
	int ans = 0;
	for (int i = 0; i < lt - lp + 1; i++) {
		for (int j = 0; j < lp; j++) {
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
	if (argc != 3) {
		printf("Usage: ./naive <PATTERN> <TEXT>\n");
		return 1;
	}
	printf("Pattern = %s\nText = %s\n", argv[1], argv[2]);
	printf("MATCHES: %d\n", naive(argv[1], argv[2]));
}
