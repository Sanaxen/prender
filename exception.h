#ifndef _EXCEPTION_H__
#define _EXCEPTION_H__

#include <string>
#include <exception>
#include <windows.h>
#include <eh.h>

inline void SETransFunc(unsigned int code, PEXCEPTION_POINTERS data)
{
	std::string desc;
	switch (code)
	{
	case EXCEPTION_ACCESS_VIOLATION:			{ desc = "EXCEPTION_ACCESS_VIOLATION";			break; }
	case EXCEPTION_ARRAY_BOUNDS_EXCEEDED:		{ desc = "EXCEPTION_ARRAY_BOUNDS_EXCEEDED";		break; }
	case EXCEPTION_BREAKPOINT:					{ desc = "EXCEPTION_BREAKPOINT";				break; }
	case EXCEPTION_DATATYPE_MISALIGNMENT:		{ desc = "EXCEPTION_DATATYPE_MISALIGNMENT";		break; }
	case EXCEPTION_FLT_DENORMAL_OPERAND:		{ desc = "EXCEPTION_FLT_DENORMAL_OPERAND";		break; }
	case EXCEPTION_FLT_DIVIDE_BY_ZERO:			{ desc = "EXCEPTION_FLT_DIVIDE_BY_ZERO";		break; }
	case EXCEPTION_FLT_INEXACT_RESULT:			{ desc = "EXCEPTION_FLT_INEXACT_RESULT";		break; }
	case EXCEPTION_FLT_INVALID_OPERATION:		{ desc = "EXCEPTION_FLT_INVALID_OPERATION";		break; }
	case EXCEPTION_FLT_OVERFLOW:				{ desc = "EXCEPTION_FLT_OVERFLOW";				break; }
	case EXCEPTION_FLT_STACK_CHECK:				{ desc = "EXCEPTION_FLT_STACK_CHECK";			break; }
	case EXCEPTION_FLT_UNDERFLOW:				{ desc = "EXCEPTION_FLT_UNDERFLOW";				break; }
	case EXCEPTION_ILLEGAL_INSTRUCTION:			{ desc = "EXCEPTION_ILLEGAL_INSTRUCTION";		break; }
	case EXCEPTION_IN_PAGE_ERROR:				{ desc = "EXCEPTION_IN_PAGE_ERROR";				break; }
	case EXCEPTION_INT_DIVIDE_BY_ZERO:			{ desc = "EXCEPTION_INT_DIVIDE_BY_ZERO";		break; }
	case EXCEPTION_INT_OVERFLOW:				{ desc = "EXCEPTION_INT_OVERFLOW";				break; }
	case EXCEPTION_INVALID_DISPOSITION:			{ desc = "EXCEPTION_INVALID_DISPOSITION";		break; }
	case EXCEPTION_NONCONTINUABLE_EXCEPTION:	{ desc = "EXCEPTION_NONCONTINUABLE_EXCEPTION";	break; }
	case EXCEPTION_PRIV_INSTRUCTION:			{ desc = "EXCEPTION_PRIV_INSTRUCTION";			break; }
	case EXCEPTION_SINGLE_STEP:					{ desc = "EXCEPTION_SINGLE_STEP";				break; }
	case EXCEPTION_STACK_OVERFLOW:				{ desc = "EXCEPTION_STACK_OVERFLOW";			break; }
	}

	printf("Structured exception:\n");
	printf("Exception code    : 0x%08x\n", code);
	printf("Exception address : 0x%08x\n", data->ExceptionRecord->ExceptionAddress);
	if (!desc.empty())
	{
		printf("Description       : %s\n", desc.c_str());
	}

	__debugbreak();

	throw std::runtime_error("<<<<<<<<<<< Aborting >>>>>>>>>>>>\n");
}


inline void set_SEfunction()
{
	_set_se_translator(SETransFunc);
}
#endif
