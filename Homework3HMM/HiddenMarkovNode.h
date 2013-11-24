#pragma once

#include <map>
#include <vector>
#include <string>

class HiddenMarkovNode
{
protected:
	std::map<char, float> emissionProbability;
	std::map<int, float> transitionProbability;
public:
	HiddenMarkovNode(void);
	~HiddenMarkovNode(void);
	void SetEmissionProbability(char letter, float probability);
	void SetTransitionProbability(int toState, float probability);
};

class HiddenMarkovMatrixNode
{
public:
	float probability;
};


class HiddenMarkovModel
{
protected:
	std::vector<HiddenMarkovNode> markovStates;
	std::vector< std::vector<HiddenMarkovMatrixNode> > matrix;
	std::string sequence;
	int beginState, endState;

public:
	HiddenMarkovModel();
	void LoadTransitionProbabilities(char * filename);
	void LoadEmissionProbabilities(char * filename);
	void SetSequence(std::string sequence);
	void SetStates(int beginState, int endState);

	float ComputeForwardAlgorithm();
};
