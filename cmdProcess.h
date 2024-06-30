#ifndef _CMDPROCESS_H
#define _CMDPROCESS_H

#include <Windows.h>

inline int cmdProcess(char* cmdline)
{
    STARTUPINFOA si;
    PROCESS_INFORMATION pi;

    ZeroMemory( &si, sizeof(si) );
    si.cb = sizeof(si);
    ZeroMemory( &pi, sizeof(pi) );


    // �G�ۂ��N������
    if( !CreateProcessA( NULL, // No module name (use command line). 
        cmdline, // Command line. 
        NULL,             // Process handle not inheritable. 
        NULL,             // Thread handle not inheritable. 
        FALSE,            // Set handle inheritance to FALSE. 
        0,                // No creation flags.
        NULL,             // Use parent's environment block.
        NULL,             // Use parent's starting directory.
        &si,              // Pointer to STARTUPINFO structure.
        &pi )             // Pointer to PROCESS_INFORMATION structure.
    )
    {
		return -1;
    }
	return 0;
}

#endif