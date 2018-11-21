#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <sstream>

using namespace std;

template<typename T>
bool comparePairs(const std::pair<T,T>& lhs, const std::pair<T,T>& rhs)
{
	return lhs.first < rhs.first;
}

typedef double Evalue;

/* an attempt is made here to match the rounding error of OrthoMCL
   by mimicking the floating point operations of mysql and perl. */
double computeScore(Evalue a, Evalue b, Evalue min)
{
	//compute the score (mysql, orthomclPairs)
	double score;
	if(a == min || b == min)
		score = (ceil(-log10(a))+ceil(-log10(b)))/2.0;
	else
		score = (log10(a) + log10(b))/(-2.0);
	
	//insert into float field (mysql, orthomclPairs)
	float s = score;

	//stringify and convert back to double (perl DBI, orthomclDumpPairsFiles)
	stringstream ss;
	ss << s;

	score = atof(ss.str().c_str());
	
	//truncate towards zero (perl, orthomclDumpPairsFiles);
	return int(score * 1000 + 0.5)/1000.0;
}

float computePML(int len_a, int len_b, vector<pair<unsigned int, unsigned int> > hsps)
{
	sort(hsps.begin(), hsps.end(), comparePairs<unsigned int>);
				
	pair<unsigned int, unsigned int> tmp_hsp = hsps.front();
			
	float pml = 0.0;	
	for(size_t i = 1; i < hsps.size(); i++)
	{
		if(hsps[i].second <= tmp_hsp.second)
			continue;
					
		if(hsps[i].first <= tmp_hsp.second)
		{
			tmp_hsp.second = hsps[i].second;
		}
		else
		{
			pml += tmp_hsp.second - tmp_hsp.first + 1;
			tmp_hsp.first = hsps[i].first;
			tmp_hsp.second = hsps[i].second;
		}
	}
				
	pml += tmp_hsp.second - tmp_hsp.first + 1;
				
	float len = len_a < len_b ? len_a : len_b;
				
	return int(pml/len * 1000 + 0.5)/10.0;
}

int main( int argc, char* argv[] )
{
	if(argc != 4)
	{
		cerr << "usage: homomcl <length-file> <e-value> <blast-file>\n\n";
		cerr << "length-file\t- tab separated table of sequence ids and their lengths\n";
		cerr << "e-value\t\t- BLAST e-value cutoff\n";
		cerr << "blast-file\t- tabular format (-outfmt 6) BLAST results file\n";
		exit(1);
	}

	ifstream lenfile(argv[1]);
	double cutoff = atof(argv[2]);
	ifstream infile(argv[3]);

	string a;	
	unsigned int len;
	map<string, size_t> indices;
	vector<unsigned int> lens;
	while(lenfile >> a >> len)
	{
		indices[a] = lens.size();
		lens.push_back(len);
	}

	string b,tmp;
	Evalue e;
	size_t qstart,qend,sstart,send;
	Evalue min = 1e-181;
	Evalue f = 0.0;
	while(infile >> a >> b >> tmp >> tmp >> tmp >> tmp >> qstart >> qend >> sstart >> send >> e >> tmp)
	{
		if(a == b)
			continue;
			
		if(e > 0.0 && e < min)
			min = e;
	}
	min /= 10;

	infile.clear();
	infile.seekg(0, ios::beg);
	
	vector<pair<unsigned int, unsigned int> > hsps;
	map<size_t, map<size_t, pair<Evalue, float> > > data;
	string c,d;
	size_t index = 0;
	float pml = 0.0;
	while(infile >> a >> b >> tmp >> tmp >> tmp >> tmp >> qstart >> qend >> sstart >> send >> e >> tmp)
	{
		if(indices.find(a) == indices.end())
		{
			cerr << "Error: could not find length for " << a << endl;
		}
		if(indices.find(b) == indices.end())
                {
                        cerr << "Error: could not find length for " << b << endl;
                }

		if(a == b)
			continue;
			
		if(a != c || b != d)
		{
			if(!c.empty())
			{
				pml = computePML(lens[indices[c]],lens[indices[d]],hsps);
				hsps.clear();
				
				if(pml >= cutoff)
				{	
					//cerr << pml << "\t" << cutoff << "\n";
					if(data.find(indices[d]) != data.end() && data[indices[d]].find(indices[c]) != data[indices[d]].end())
					{
						pair<Evalue,float> pt = data[indices[d]][indices[c]];
						
						double score = computeScore(pt.first, f, min);
		
						cout << c << "\t" << d << "\t" << score << "\n";
	
						data[indices[d]].erase(indices[c]);
					}
					else
					{
						data[indices[c]][indices[d]] = pair<Evalue,float>(f,pml);
					}
				}
			}
			
			if(indices.find(a) == indices.end())
				indices[a] = index++;
		
			if(indices.find(b) == indices.end())
				indices[b] = index++;
		
			c = a;
			d = b;
		}
		
		pair<unsigned int, unsigned int> hsp, hspa, hspb;
		
		hspa.first = qstart;
		hspa.second = qend;
		
		hspb.first = sstart;
		hspb.second = send;
		
		hsp = lens[indices[a]] < lens[indices[b]] ? hspa : hspb;
	
		hsps.push_back(hsp);
		
		if(hsps.size() == 1)
		{
			f = min;
			if(e > 0.0)
				f = e;
		}
	}
	
	/* don't forget the last one */

	pml = computePML(lens[indices[c]],lens[indices[d]],hsps);
	
	if(pml >= cutoff)
	{
		if(data.find(indices[d]) != data.end() && data[indices[d]].find(indices[c]) != data[indices[d]].end())
		{
			pair<Evalue,float> pt = data[indices[d]][indices[c]];
			
			double score = computeScore(pt.first, f, min);
		
			cout << c << "\t" << d << "\t" << score << "\n";
	
			data[indices[d]].erase(indices[c]);
		}
		else
		{
			data[indices[c]][indices[d]] = pair<Evalue,float>(f,pml);
		}
	}

	vector<string> seqs(indices.size(), "");
	for(map<string, size_t>::iterator it = indices.begin(); it != indices.end();)
	{
		seqs[it->second] = it->first;
		indices.erase(it++);
	}
	/* find orphan edges */
	/*cout << "-------------------------------------------------------------------------\n";
	for(map<size_t, map<size_t, pair<Evalue,float> > >::iterator it = data.begin(); it != data.end(); it++)
	{
		for(map<size_t, pair<Evalue,float> >::iterator jt = it->second.begin(); jt != it->second.end(); jt++)
		{
			size_t i = it->first;
			size_t j = jt->first;

			pair<Evalue,float> pt = data[i][j];
			
			double score;
			if(pt.first == min)
				score = ceil(-log10(pt.first));
			else
				score = -log10(pt.first);
				
			float s = score;        		
			
			stringstream ss;
			ss << s;
			
			score = atof(ss.str().c_str());
			score = int(score * 1000.0 + .5)/1000.0;
			
			cout << seqs[i] << "\t" << seqs[j] << "\t" << score << "\n";
		}
	}*/
}
