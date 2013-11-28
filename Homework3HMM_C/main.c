#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

//#define FORWARD_ALGORITHM
//#define VITERBI_ALGORITHM
#define BACKWARD_ALGORITHM

#ifdef FORWARD_ALGORITHM
char * usage = "usage: ./doForward transitions.txt emissions.txt begin-state end-state sequence \n";
#endif

#ifdef VITERBI_ALGORITHM
char * usage = "usage: ./doBackward transitions.txt emissions.txt begin-state end-state sequence \n";
#endif

#ifdef BACKWARD_ALGORITHM
char * usage = "usage: ./doViterbi transitions.txt emissions.txt begin-state end-state sequence \n";
#endif

//A,C,G,T   (26 for proteins)
#define alphabetSize 4
#define maxStates 6


//globals
char * alphabetFormat = "ACGT";
int numStates = 0;
int beginState = 0, endState = 0;
char * sequence = NULL;
char * convertedSequence = NULL;
size_t sequenceLength = 0;

float transitionProbabilities[maxStates-1][maxStates];	//-1 because end state does not transition
float emissionProbabilites[maxStates-1][alphabetSize];	//-1 because end state does not emit
float * probabilityMatrix = NULL;	//primary data structure (width = states, height = sequence)
int * viterbiBackPointers = NULL;	//Viterbi algorithm data structure to keep track of back pointers

int FindChar(char * string, char c, int length)
{
	char * found = (char*)memchr(string, c, length);
	if(found)
		return found - string;
	return -1;
}

//Find the forward probability of an individual cell
//shouldEmit specifies whether this is a normal node(true), or an end state -> don't emit
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

//Compute an individual cell's Viterbi probability and backpointer
float ViterbiProbability(int x, int y, int shouldEmit)
{
	int toState = x;
	float maxProbability = 0.0f;
	int fromState, maxState = x;
	float transitionProbability, emissionProbability, previousProbability, stateProbability;

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
		stateProbability = transitionProbability * emissionProbability * previousProbability;
		if(stateProbability > maxProbability)
		{
			maxProbability = stateProbability;
			maxState = fromState;
		}
	}
	probabilityMatrix[y * numStates + x] = maxProbability;
	viterbiBackPointers[y * numStates + x] = maxState;
	return maxProbability;
}


//Display the Viterbi Path
void ViterbiPath(int endX, int endY)
{
	int y, x;	//y = sequence, x = state
	size_t statesSize;
	int * states;
	
	statesSize = (sequenceLength + 2) * sizeof(int);
	states = (int*)malloc(statesSize);

	x = endX;
	for(y = endY; y > 0; y--)
	{
		x = viterbiBackPointers[y * numStates + x];
		states[y] = x;
	}
	
	printf("Viterbi path: \n");
	printf("0 \n");
	for(y = 2; y <= endY; y++)
	{
		printf("%d -> %c \n", states[y], sequence[y - 2]);
	}


	free(states);
}

//Compute a cell's backward probability
//shouldEmit specifies whether this is a normal node(true), or an end state -> don't emit
float BackwardProbability(int x, int y, int shouldEmit)
{
	int fromState = x;
	float probability = 0.0f;
	int toState;
	float transitionProbability, emissionProbability, previousProbability;

	for(toState = 0; toState < numStates - 1; toState++)
	{
		transitionProbability = transitionProbabilities[fromState][toState];
		emissionProbability = 1.0f;
		previousProbability = 1.0f;
		if(shouldEmit)
		{
			char letter = convertedSequence[y];
			assert(letter < alphabetSize);
			emissionProbability = emissionProbabilites[toState][letter];
			previousProbability = probabilityMatrix[ (y+1) * numStates + toState];
		}
		probability += transitionProbability * emissionProbability * previousProbability;
	}
	probabilityMatrix[y * numStates + x] = probability;
	return probability;
}


//Computes the Forward Algorithm Matrix as well as Viterbi
int DoForward(int mode)
{
	int x, y, inverseX, inverseY, totalSequenceLength;
	size_t probabilityMatrixSize, viterbiBackPointersSize;

	probabilityMatrixSize = (numStates) * (sequenceLength + 2) * sizeof(float);
	probabilityMatrix = (float*)malloc(probabilityMatrixSize);
	memset(probabilityMatrix, 0, probabilityMatrixSize);

	//create viterbi data structure
	if(mode == 1)
	{
		viterbiBackPointersSize = (numStates) * (sequenceLength + 2) * sizeof(int);
		viterbiBackPointers = (int*)malloc(viterbiBackPointersSize);
		memset(viterbiBackPointers, 0, viterbiBackPointersSize);
	}


	//Fill in edge cases
	switch(mode)
	{
	case 0:	//forward
	case 1:	//viterbi
		probabilityMatrix[0] = 1.0f;
		for(x = 1; x < numStates; x++)
			probabilityMatrix[x] = 0.0f;
		for(y = 1; y < sequenceLength + 1; y++)
			probabilityMatrix[y * numStates] = 0.0f;
		break;
	case 2:	//backwards
		totalSequenceLength = (sequenceLength+2);
		probabilityMatrix[totalSequenceLength * numStates - 1] = 1.0f;
		totalSequenceLength--;
		for(x = 0; x < numStates - 1; x++)
			probabilityMatrix[(totalSequenceLength-1) * numStates + x] = transitionProbabilities[x][numStates-1];
		for(y = 0; y < totalSequenceLength; y++)
			probabilityMatrix[y * numStates + numStates - 1] = 0.0f;
		//handle no-emit case
		//printf("Beta for state %d time %d: %f \n", x, y-1, BackwardProbability(numStates- 1, sequenceLength + 1, 0) );
		break;
	default:
		assert(0);
	}


	//Fill in the Matrix
	for(y = 1; y < sequenceLength + 1; y++)
	{
		for(x = 1; x < numStates - 1; x++)
		{
			switch(mode)
			{
			case 0:	//Forward Algorithm
				printf("Alpha for state %d time %d: %f \n", x, y-1, ForwardProbability(x, y, 1) );
				break;
			case 1:	//Viterbi Algorithm
				ViterbiProbability(x, y, 1);
				printf("Viterbi for state %d time %d: %f maxstate %d \n", x, y-1, 
					probabilityMatrix[y * numStates + x], viterbiBackPointers[y * numStates + x] );
				break;
			case 2:	//Backward Algorithm
				if(y == 1) continue;	//skip the first row, handled in edge case section
				inverseX = numStates - x - 1;
				inverseY = sequenceLength + 1 - y;
				printf("Beta for state %d time %d: %f \n", x, y-1, BackwardProbability(inverseX, inverseY, 1) );
				break;
			default:
				assert(0);
			}
		}
	}
	
	switch(mode)
	{
	case 0:	//Forward Algorithm
		printf("Forward probability: %f \n", ForwardProbability(numStates - 1, sequenceLength + 1, 0) );
		break;
	case 1:	//Viterbi Algorithm
		printf("Viterbi probability: %f \n", ViterbiProbability(numStates - 1, sequenceLength + 1, 0) );
		ViterbiPath(numStates - 1, sequenceLength + 1);
		break;
	case 2:	//Backward Algorithm
		printf("Backward probability: %f \n", BackwardProbability(0, 0, 1) );
		break;
	default:
		assert(0);
	}
	
	//Cleanup
	free(probabilityMatrix);
	if(mode == 1)
		free(viterbiBackPointers);

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
	convertedSequence = (char*)malloc(sequenceLength + 2);
	for(i = 0; i < sequenceLength; i++)
	{
		letter = sequence[i];
		sequenceIndex = FindChar(alphabetFormat, letter, alphabetSize);
		assert(sequenceIndex >= 0);		//letter not in alphabet
		convertedSequence[i] = sequenceIndex;
	}
	convertedSequence[sequenceLength] = 0;
}



int main(int argc, char * argv[])
{
	if(argc < 6) {
		printf(usage);
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
	DoForward(0);
#endif
#ifdef VITERBI_ALGORITHM
	DoForward(1);
#endif
#ifdef BACKWARD_ALGORITHM
	DoForward(2);
#endif

	//Cleanup

	free(convertedSequence);

	return 0;
}
