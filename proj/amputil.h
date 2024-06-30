#ifndef _AMPUTIL_H_

#define _AMPUTIL_H_
#include <amp.h>
#include <amp_math.h>
using namespace concurrency;

inline std::string WStringToString
(
std::wstring oWString
)
{
	// wstring → SJIS
	int iBufferSize = WideCharToMultiByte(CP_OEMCP, 0, oWString.c_str()
		, -1, (char *)NULL, 0, NULL, NULL);

	// バッファの取得
	CHAR* cpMultiByte = new CHAR[iBufferSize];

	// wstring → SJIS
	WideCharToMultiByte(CP_OEMCP, 0, oWString.c_str(), -1, cpMultiByte
		, iBufferSize, NULL, NULL);

	// stringの生成
	std::string oRet(cpMultiByte, cpMultiByte + iBufferSize - 1);

	// バッファの破棄
	delete[] cpMultiByte;

	// 変換結果を返す
	return(oRet);
}

inline accelerator acceleratorInfo()
{
	std::vector<accelerator> accs = accelerator::get_all();
	for (int i = 0; i < accs.size(); i++)
	{
		printf("[%d]description:%s\n", i, WStringToString(accs[i].description).c_str());
		printf("[%d]device_path:%s\n", i, WStringToString(accs[i].device_path).c_str());
		printf("[%d]dedicated_memory:%d\n", i, accs[i].dedicated_memory);
		printf("[%d]supports_double_precision:%s\n\n",
			i, accs[i].supports_double_precision ? "true" : "false");
	}

	accelerator default_acc;
	printf("default_acc.description:%s\n", WStringToString(default_acc.description).c_str());
	printf("default_acc.device_path:%s\n", WStringToString(default_acc.device_path).c_str());
	printf("default_acc.dedicated_memory:%d\n", default_acc.dedicated_memory);
	printf("default_acc.supports_double_precision:%s\n\n",
		default_acc.supports_double_precision ? "true" : "false");

	return default_acc;
}

#endif