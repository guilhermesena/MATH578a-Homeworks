#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>

using namespace std;

class Sequence {
public:
	string name;
	string seq;

	Sequence() {
		name = "";
		seq = "";
	}

	Sequence(string arg_name) {
		name = arg_name;
		seq = "";
	}
};

class Alignment {
private:
	static const int SCORE_MATCH = 3;
	static const int SCORE_MISMATCH = -1;
	static const int SCORE_GAP = -3;

	Sequence sa, sb;
	size_t la, lb;

	int ** score;
	pair<int, int> ** prev;

public:
	Alignment(Sequence arg_sa, Sequence arg_sb) {
		sa = arg_sa;
		sb = arg_sb;

		la = sa.seq.size() + 1;
		lb = sb.seq.size() + 1;

		score = new int*[la + 1];
		prev = new pair<int, int>*[la];
		for (int i = 0; i < la; i++) {
			score[i] = new int[lb + 1];
			prev[i] = new pair<int, int> [lb + 1];
		}
	}

	void Needleman_Wunsch() {
		for (int i = 0; i < la; i++) {
			score[i][0] = i * SCORE_GAP;
			prev[i][0] = make_pair(i, 0);
		}

		for (int i = 0; i < lb; i++) {
			score[0][i] = i * SCORE_GAP;
			prev[0][i] = make_pair(0, i);
		}

		for (int i = 1; i < la; i++) {
			for (int j = 1; j < lb; j++) {
				int t1 = score[i - 1][j] + SCORE_GAP;
				int t2 = score[i][j - 1] + SCORE_GAP;
				int t3 = score[i - 1][j - 1]
						+ ((sa.seq[i] == sb.seq[j]) ?
								SCORE_MATCH : SCORE_MISMATCH);

				score[i][j] = t1;
				prev[i][j] = make_pair(i - 1, j);
				if (t2 > score[i][j]) {
					score[i][j] = t2;
					prev[i][j] = make_pair(i, j - 1);
				}
				if (t3 > score[i][j]) {
					score[i][j] = t3;
					prev[i][j] = make_pair(i - 1, j - 1);
				}
			}
		}
	}

	void print_alignment() {
		pair<int, int> cur = make_pair(la - 1, lb - 1);
		vector<char> a_backwards, b_backwards;

		while (cur != prev[cur.first][cur.second]) {
			if (cur.first == prev[cur.first][cur.second].first) {
				a_backwards.push_back(sa.seq[cur.first - 1]);
				b_backwards.push_back('-');
			}

			else if (cur.second == prev[cur.first][cur.second].second) {
				a_backwards.push_back('-');
				b_backwards.push_back(sa.seq[cur.second - 1]);
			} else {
				a_backwards.push_back(sa.seq[cur.first - 1]);
				b_backwards.push_back(sb.seq[cur.second - 1]);
			}
			cur = prev[cur.first][cur.second];
		}

		reverse(a_backwards.begin(), a_backwards.end());
		printf("%s\n", string(a_backwards.begin(), a_backwards.end()).c_str());

		reverse(b_backwards.begin(), b_backwards.end());
		printf("%s\n", string(b_backwards.begin(), b_backwards.end()).c_str());

		printf("Score = %d\n\n", score[la - 1][lb - 1]);
	}
};

class Aligner {
public:

	void read(ifstream &inp) {
		string word_read;
		vector<Sequence> v;
		int nseq = -1;
		while (!inp.eof()) {
			inp >> word_read;

			if (word_read[0] == '>') {
				++nseq;
				string title = word_read;
				getline(inp, word_read);
				title.append(word_read);
				v.push_back(Sequence(title));
			} else if (word_read.size() > 0) {
				v[nseq].seq.append(word_read);
			}
		}
		inp.close();

		int len = v.size();
		for (int i = 0; i < len; i++) {
			printf("Sequence %d: %s\n%s\n\n", i + 1, v[i].name.c_str(),
					v[i].seq.c_str());
		}

		for (int i = 0; i < len; i++) {
			for (int j = 0; j < len; j++) {

				printf("Aligning sequences %d and %d...\n", i + 1, j + 1);

				Alignment a(v[i], v[j]);
				a.Needleman_Wunsch();
				a.print_alignment();
			}
		}

	}

	int run(int argc, char** argv) {
		if (argc != 2) {
			printf("Please specify exactly 1 input file\n");
			return 1;
		}

		ifstream inp;
		inp.open(argv[1]);

		if (!inp.is_open()) {
			printf("Invalid input file\n");
			return 1;
		}

		read(inp);
		return 0;
	}
};

int main(int argc, char** argv) {
	Aligner a;
	a.run(argc, argv);
}
