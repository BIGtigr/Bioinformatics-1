#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define FORWARD_ALGORITHM
//#define VITERBI_ALGORITHM
//#define BACKWARD_ALGORITHM

//A,C,G,T   (26 for proteins)
#define alphabetSize 4
#define maxStates 6


//globals
char * alphabetFormat = "ACGT";
int numStates = 0;
int beginState = 0, endState = 0;
char * sequence;
char * convertedSequence;
size_t sequenceLength = 0;

float transitionProbabilities[maxStates-1][maxStates];	//-1 because end state does not transition
float emissionProbabilites[maxStates-1][alphabetSize];	//-1 because end state does not emit
float * probabilityMatrix;	//primary data structure (width = states, height = sequence)

int FindChar(char * string, char c, int length)
{
	char * found = (char*)memchr(string, c, length);
	if(found)
		return found - string;
	return -1;
}


float ForwardProbability(int x, int y, int shouldEmit)
{
	int toState = x;
	float probability = 0.0f;
	int fromState;
	float transitionProbability, emissionProbability, previousProbability;

	for(fromState = 0; fromState < numStates - 1; fromState++)
	{
		transitionProbability = transitionProbabilities[fromState][toState];
		emissionProbability = 1.0f;
		if(shouldEmit)
		{
			char letter = convertedSequence[y - 1];
			assert(letter < alphabetSize);
			emissionProbability = emissionProbabilites[toState][letter];
		}
		previousProbability = probabilityMatrix[ (y - 1) * numStates + fromState];
		probability += transitionProbability * emissionProbability * previousProbability;
	}
	probabilityMatrix[y * numStates + x] = probability;
	return probability;
}


int DoForward()
{
	int x, y;

	size_t probabilityMatrixSize = (numStates) * (sequenceLength + 1) * sizeof(float);
	probabilityMatrix = (float*)malloc(probabilityMatrixSize);
	memset(probabilityMatrix, 0, probabilityMatrixSize);

	//Fill in edge cases
	probabilityMatrix[0] = 1.0f;
	for(x = 1; x < numStates; x++)
		probabilityMatrix[x] = 0.0f;
	for(y = 1; y < sequenceLength + 1; y++)
		probabilityMatrix[y * numStates] = 0.0f;

	//Fill in the Matrix
	for(y = 1; y < sequenceLength + 1; y++)
	{
		for(x = 1; x < numStates - 1; x++)
		{
			printf("Alpha for state %d time %d: %f \n", x, y, ForwardProbability(x, y, 1) );
		}
	}
	
	printf("Forward probability: %f \n", ForwardProbability(numStates - 1, sequenceLength + 1, 0) );

	free(probabilityMatrix);

	return 0;
}


float ViterbiProbability(int x, int y, int shouldEmit)
{
	int toState = x;
	float probability = 0.0f;
	int fromState;
	float transitionProbability, emissionProbability, previousProbability;

	for(fromState = 0; fromState < numStates - 1; fromState++)
	{
		transitionProbability = transitionProbabilities[fromState][toState];
		emissionProbability = 1.0f;
		if(shouldEmit)
		{
			char letter = convertedSequence[y - 1];
			assert(letter < alphabetSize);
			emissionProbability = emissionProbabilites[toState][letter];
		}
		previousProbability = probabilityMatrix[ (y - 1) * numStates + fromState];
		probability += transitionProbability * emissionProbability * previousProbability;
	}
	probabilityMatrix[y * numStates + x] = probability;
	return probability;
}


int DoViterbi()
{
	int x, y;

	size_t probabilityMatrixSize = (numStates) * (sequenceLength + 1) * sizeof(float);
	probabilityMatrix = (float*)malloc(probabilityMatrixSize);
	memset(probabilityMatrix, 0, probabilityMatrixSize);

	//Fill in edge cases
	probabilityMatrix[0] = 1.0f;
	for(x = 1; x < numStates; x++)
		probabilityMatrix[x] = 0.0f;
	for(y = 1; y < sequenceLength + 1; y++)
		probabilityMatrix[y * numStates] = 0.0f;

	//Fill in the Matrix
	for(y = 1; y < sequenceLength + 1; y++)
	{
		for(x = 1; x < numStates - 1; x++)
		{
			printf("Alpha for state %d time %d: %f \n", x, y, ForwardProbability(x, y, 1) );
		}
	}
	
	printf("Forward probability: %f \n", ForwardProbability(numStates - 1, sequenceLength + 1, 0) );

	free(probabilityMatrix);

	return 0;
}

void ReadTransitionProbabilities(char * filename)
{
	FILE *file;
	int fromState, toState;
	float probability;

	file = fopen(filename, "r");
	assert(file);

	while( fscanf(file, "%d %d %f", &fromState, &toState, &probability) != EOF )
		transitionProbabilities[fromState][toState] = probability;

	fclose(file);
}

void ReadEmissionProbabilities(char * filename)
{
	FILE *file;
	int state, index;
	char letter;
	float probability;

	file = fopen(filename, "r");
	assert(file);
	
	while( fscanf(file, "%d %c %f", &state, &letter, &probability) != EOF ) {
		assert(state < maxStates);

		index = FindChar(alphabetFormat, letter, alphabetSize);
		assert(index >= 0);		//letter not in alphabet
		emissionProbabilites[state][index] = probability;

		if(numStates < state)	//update the actual number of states
			numStates = state;
	}

	numStates += 2;		//beginning and end states (not in emission)

	fclose(file);
}

//Convert {'A','C','G','T'} to {0, 1, 3, 4}
//Can be easily modified for proteins by changing alphabetFormat, alphabetSize
void ConvertSequence()
{
	int i, sequenceIndex;
	char letter;

	assert(sequence && sequenceLength > 0);
	convertedSequence = (char*)malloc(sequenceLength + 1);
	for(i = 0; i < sequenceLength; i++)
	{
		letter = sequence[i];
		sequenceIndex = FindChar(alphabetFormat, letter, alphabetSize);
		assert(sequenceIndex >= 0);		//letter not in alphabet
		convertedSequence[i] = sequenceIndex;
	}
	convertedSequence[sequenceLength] = 0;
}



char * usageForward = "usage: ./doForward transitions.txt emissions.txt begin-state end-state sequence \n";
char * usageBackward = "usage: ./doBackward transitions.txt emissions.txt begin-state end-state sequence \n";
char * usageViterbi = "usage: ./doViterbi transitions.txt emissions.txt begin-state end-state sequence \n";
int main(int argc, char * argv[])
{
	if(argc < 6) {
		#ifdef FORWARD_ALGORITHM
			printf(usageForward);
		#endif
		#ifdef VITERBI_ALGORITHM
			printf(usageViterbi);
		#endif
		#ifdef BACKWARD_ALGORITHM
			printf(usageBackward);
		#endif
		return 1;
	}

	//Read in Transition Probabilities
	ReadTransitionProbabilities(argv[1]);
	//Read in Emission Probabilities
	ReadEmissionProbabilities(argv[2]);

	sscanf(argv[3], "%d", &beginState);
	sscanf(argv[4], "%d", &endState);
	sequence = argv[5];
	sequenceLength = strlen(sequence);

	ConvertSequence();

#ifdef FORWARD_ALGORITHM
	DoForward();
#endif
#ifdef VITERBI_ALGORITHM
	DoViterbi();
#endif
#ifdef BACKWARD_ALGORITHM
	DoBackward();
#endif

	//Cleanup
	//free(convertedSequence);
	
	return 0;
}
