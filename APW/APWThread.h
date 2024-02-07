#pragma once

#include <thread>
#include <vector>
#include <atomic>

#include "Options.h"

class APWFrame;

class APWThread
{
public:
	APWThread(const Options& options, APWFrame* frame);
	~APWThread();


	void Start();
	void join();

	void Terminate();


	const Options& m_options;

	APWFrame* m_frame;

	std::vector<std::vector<double>> results;

private:
	void Calculate();

	std::thread mThread;
	std::atomic_bool terminate;
};

