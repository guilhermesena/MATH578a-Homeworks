#include <iostream>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <algorithm>
#include <string.h>
#include <map>
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
	static const int INF = 1 << 30;
	Sequence sa, sb;
	map<pair<int,int> , pair<int,int> > prev;
	int la, lb;

	int char_dist (char a, char b) {
		return (a==b) ? DIST_MATCH: DIST_MISMATCH;
	}

public:
	Alignment(Sequence arg_sa, Sequence arg_sb) {
		//Guarantees sa always has the smaller length
		if(arg_sa.seq.size() < arg_sb.seq.size()) {
			sa = arg_sa;
			sb = arg_sb;
		}
		else {
			sa = arg_sb;
			sb = arg_sa;
		}

		//la < lb
		la = sa.seq.size() + 1;
		lb = sb.seq.size() + 1;

	}

	int Banded_DP (int k) {
		//la rows, (lb - la +1) + 2k columns
		int delta = lb - la + 1 + 2*k; 
		int **dist = new int*[la];

		for(int i = 0; i < la; i++) {
			dist[i] = new int[delta];
		}
		prev.clear();
	
		for(int i = 0; i < la; i++) {
			int left = max(0, i - k);
			int right = min(lb, i - k + delta);

			for(int j = 0; j < right - left; j++) {
				int u = i, v = j + left; 
				if(u == 0) {
					dist[i][j] = v*DIST_GAP;
					prev[make_pair(u,v)] = make_pair(u,v);
					continue;
				}	

				if(v == 0) {
					dist[i][j] = u*DIST_GAP;
					prev[make_pair(u,v)] = make_pair(u,v);
					continue;
				}

				int t1=INF,t2=INF,t3;

				if(j>0) 
					t1 = dist[i][j-1] + DIST_GAP;
				if(j < delta-1) 
					t2 = dist[i-1][j+((left == 0)?0:1)] + DIST_GAP;
				
				t3 = dist[i-1][j - ((left == 0)?1:0)] + char_dist(sa.seq[u-1], sb.seq[v-1]);

				dist[i][j] = t3;
				prev[make_pair(u,v)] = make_pair(u-1,v-1);
				
				if(t1 < dist[i][j]) {
					dist[i][j] = t1;
					prev[make_pair(u,v)] = make_pair(u,v-1);
				}

				if(t2 < dist[i][j]) {
					dist[i][j] = t2;
					prev[make_pair(u,v)] = make_pair(u-1,v);
				}
			
			}
		}
		return dist[la-1][min(lb-1,delta-k-1)];
	}

	void print_alignment(int dist) {
		pair<int, int> cur = make_pair(la - 1, lb - 1);
		vector<char> a_backwards, b_backwards;

		while (cur != prev[cur]) {
			if (cur.first == prev[cur].first) {
				a_backwards.push_back(sa.seq[cur.first - 1]);
				b_backwards.push_back('-');
			}

			else if (cur.second == prev[cur].second) {
				a_backwards.push_back('-');
				b_backwards.push_back(sa.seq[cur.second - 1]);
			} else {
				a_backwards.push_back(sa.seq[cur.first - 1]);
				b_backwards.push_back(sb.seq[cur.second - 1]);
			}
			cur = prev[cur];
		}

		reverse(a_backwards.begin(), a_backwards.end());
		printf("%s\n", string(a_backwards.begin(), a_backwards.end()).c_str());

		reverse(b_backwards.begin(), b_backwards.end());
		printf("%s\n", string(b_backwards.begin(), b_backwards.end()).c_str());

		printf("Distance = %d\n\n", dist);
	}
};

class Aligner {
public:
	void read(ifstream &inp, vector<Sequence> &v) {
		string word_read;
		int nseq = -1;

		inp >> word_read;
		while (!inp.eof()) {
			if (word_read[0] == '>') {
				++nseq;
				string title = word_read;
				getline(inp, word_read);
				title.append(word_read);
				v.push_back(Sequence(title));
			} else if (word_read.size() > 0) {
				v[nseq].seq.append(word_read);
			}
			inp >> word_read;
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
				printf("Aligning sequences %d and %d...\n",i+1,j+1);
				Alignment a(v[i], v[j]);
				int dist = a.Banded_DP(1),k=1;
				
				while(dist > k) {
					k <<= 1;
					dist = a.Banded_DP(k);
				}
				a.print_alignment(dist);
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
