#include "HiddenMarkovNode.h"

#include <fstream>
#include <assert.h>

using namespace std;

// HIDDEN MARKOV NODE =========================================================
HiddenMarkovNode::HiddenMarkovNode(void)
{
}


HiddenMarkovNode::~HiddenMarkovNode(void)
{
}

void HiddenMarkovNode::SetEmissionProbability(char letter, float probability)
{
	emissionProbability[letter] = probability;
}

void HiddenMarkovNode::SetTransitionProbability(int toState, float probability)
{
	transitionProbability[toState] = probability;
}

// HIDDEN MARKOV MODEL ========================================================
HiddenMarkovModel::HiddenMarkovModel()
{
}

void HiddenMarkovModel::LoadTransitionProbabilities(char * filename)
{
	ifstream file(filename);

	while(file.good()) {
		unsigned int fromState, toState;
		float probability;

		file >> fromState;
		file >> toState;
		file >> probability;

		assert(fromState >= 0);
		assert(toState > 0);
		assert(probability >= 0.0f && probability <= 1.0f);

		if(markovStates.size() < fromState + 1)
			markovStates.resize(fromState + 1);
		/*if(markovStates.size() < toState + 1)
			markovStates.resize(toState + 1);*/

		markovStates[fromState].SetTransitionProbability(toState, probability);
	}
}

void HiddenMarkovModel::LoadEmissionProbabilities(char * filename)
{
	ifstream file(filename);

	while(file.good()) {
		unsigned int state;
		char letter;
		float probability;

		file >> state;
		file >> letter;
		file >> probability;

		assert(state > 0);
		assert(probability >= 0.0f && probability <= 1.0f);
		assert(letter);

		if(markovStates.size() < state + 1)
			markovStates.resize(state + 1);

		markovStates[state].SetEmissionProbability(letter, probability);
	}

}


void HiddenMarkovModel::SetSequence(string sequence)
{
	this->sequence = sequence;
}


void HiddenMarkovModel::SetStates(int beginState, int endState)
{
	this->beginState = beginState;
	this->endState = endState;
}

float HiddenMarkovModel::ComputeForwardAlgorithm()
{
	return 0.0f;
}
