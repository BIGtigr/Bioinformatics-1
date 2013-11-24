
#include "HiddenMarkovNode.h"

#include <iostream>

using namespace std;

char * usage = "usage: ./doForward transitions.txt emissions.txt begin-state end-state sequence";

void ReadTransitions(char * filename)
{

}


int main(int argc, char * argv[])
{

	if(argc < 6) {
		cout << usage;
		return 1;
	}

	char * transitionFilename = argv[1];
	char * emissionFilename = argv[2];
	int beginState = atoi(argv[3]);
	int endState = atoi(argv[4]);
	char * sequence = argv[5];

	HiddenMarkovModel hmm;
	hmm.LoadEmissionProbabilities(emissionFilename);
	hmm.LoadTransitionProbabilities(transitionFilename);
	hmm.SetSequence(sequence);
	hmm.SetStates(beginState, endState);

	return 0;
}
