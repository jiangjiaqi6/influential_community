/*
 * Timer.h
 *
 *  Created on: 5Dec.,2017
 *      Author: Lijun Chang
 *      Email: ljchang@outlook.com
 */

#ifndef TIMER_H_
#define TIMER_H_

#include "Defines.h"


//for linux
// class Timer {
// public:
// 	Timer() { m_start = timestamp(); }
// 	void restart() { m_start = timestamp(); }
// 	long long elapsed() { return timestamp() - m_start; }

// private:
// 	long long m_start;

// 	// Returns a timestamp ('now') in microseconds
// 	long long timestamp() {
// 		struct timeval tp;
// 		gettimeofday(&tp, nullptr);
// 		return ((long long)(tp.tv_sec))*1000000 + tp.tv_usec;
// 	}
// };

//for windows
class Timer {
public:
	Timer() { m_start = clock(); }
	void restart() { m_start = clock(); }
	long long elapsed() { return (clock() - m_start) / CLOCKS_PER_SEC; }

private:
	long long m_start;
};

#endif /* TIMER_H_ */