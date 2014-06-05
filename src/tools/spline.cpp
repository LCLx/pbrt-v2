#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <core/splines.h>

using namespace std;

static void usage() {
    fprintf(stderr, "usage: spline <input> <output.txt>\n");
    exit(1);
}

int main(int argc, char *argv[]) 
{
    const char *outfile = NULL;
    const char *infile = NULL;

    if (argc < 3) usage();
    infile = argv[1];
    outfile = argv[2];

    Catmull_Rom spline;
    ifstream input(infile);
    while(!input.eof()){
        float x,y;
        input>>y>>x;
        spline.Add_Sample(x,y);
    }
    input.close();
    spline.Build();
    ofstream output(outfile);
    float start_x=spline.samples[0].first;
    float step=(spline.samples[spline.samples.size()-1].first-spline.samples[0].first)/(float)(spline.samples.size()-1);
    for(size_t i=0;i<spline.samples.size();i++){
        float x=start_x+step*(float)i;
        output<<x<<" "<<spline.Evaluate(x)<<endl;
    }
    output.close();

    return 0;
}
