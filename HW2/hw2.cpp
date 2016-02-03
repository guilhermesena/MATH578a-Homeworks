#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <string.h>

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
	static const int DIST_MATCH = 0;
	static const int DIST_MISMATCH = 1;
	static const int DIST_GAP = 1;
	static const int INF = 2000000;
	Sequence sa, sb;
	int la, lb;

	int ** dist;
	pair<int, int> ** prev;

public:
	Alignment(Sequence arg_sa, Sequence arg_sb) {
		sa = arg_sa;
		sb = arg_sb;

		la = sa.seq.size() + 1;
		lb = sb.seq.size() + 1;

		dist = new int*[la + 1];
		prev = new pair<int, int>*[la];
		for (int i = 0; i < la; i++) {
			dist[i] = new int[lb + 1];
			prev[i] = new pair<int, int> [lb + 1];
		}
	}

	int Banded_DP (int k) {
		int delta = k;
		for (int i = 0; i <= delta; i++) {
			dist[i][0] = i * DIST_GAP;
			prev[i][0] = make_pair(i, 0);
		}

		for (int i = 0; i <= delta; i++) {
			dist[0][i] = i * DIST_GAP;
			prev[0][i] = make_pair(0, i);
		}

		for (int i = 1; i < la; i++) {
			int left = max(1, i - delta);
			int right = min(lb-1, i+delta);
			for (int j = left; j <= right; j++) {
				int t1 = dist[i - 1][j] + DIST_GAP;
				int t2 = dist[i][j - 1] + DIST_GAP;
				int t3 = dist[i - 1][j - 1]
						+ ((sa.seq[i-1] == sb.seq[j-1]) ?
								DIST_MATCH : DIST_MISMATCH);
				if(abs(i-1-j)>2*delta) t1 = INF;
				if(abs(i-j+1)>2*delta) t2 = INF;
				if(abs(i-j) > 2*delta) t3 = INF;

				dist[i][j] = t3;
				prev[i][j] = make_pair(i - 1, j-1);
				if (t2 < dist[i][j]) {
					dist[i][j] = t2;
					prev[i][j] = make_pair(i, j - 1);
				}
				if (t1 < dist[i][j]) {
					dist[i][j] = t1;
					prev[i][j] = make_pair(i - 1, j);
				}
			}
		}
		return dist[la-1][lb-1];
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

		printf("Distance = %d\n\n", dist[la - 1][lb - 1]);
	}
};

class Aligner {
public:
	void read(ifstream &inp, vector<Sequence> &v) {
		string word_read;
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
	}

	void print_sequences(vector<Sequence> v) {
		int len = v.size();
		for (int i = 0; i < len; i++) {
			printf("Sequence %d: %s\n%s\n\n", i + 1, v[i].name.c_str(),
					v[i].seq.c_str());
		}
	}

	void run_all_pairwise_alignments(vector<Sequence> v) {
		int len = v.size();
		for (int i = 0; i < len; i++) {
			for (int j = 0; j < len; j++) {
				printf("Aligning sequences %d and %d...\n", i + 1, j + 1);

				Alignment a(v[i], v[j]);
				int opt;
				for(int k = 1; a.Banded_DP(k) > k; k <<= 1);
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

		vector<Sequence> v;
		read(inp, v);
		print_sequences(v);
		run_all_pairwise_alignments(v);

		return 0;
	}
};

int main(int argc, char** argv) {
	Aligner a;
	return a.run(argc, argv);
}
