#include "APWThread.h"

#include "APWFrame.h"

APWThread::APWThread(const Options& options, APWFrame* frame)
	: m_options(options), m_frame(frame), terminate(false)
{
}


APWThread::~APWThread()
{
	join();
}

void APWThread::Start()
{
	mThread = std::thread([this]() {
		Calculate();
		});
}

void APWThread::join()
{
	if (mThread.joinable()) mThread.join();
}


void APWThread::Calculate()
{
	if (0 == m_options.method)
		results = m_frame->bandStructureAPW.Compute(terminate, m_options);
	else
		results = m_frame->bandStructureLAPW.Compute(terminate, m_options);

	--m_frame->runningThreads;
}

void APWThread::Terminate()
{
	terminate = true;
}
