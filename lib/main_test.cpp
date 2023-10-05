#include "test\main_test.h"

static bool debugEnabled = false;

bool isDebugEnabled() 
{
    return debugEnabled;
}

void setDebugEnabled(bool enable) 
{
    debugEnabled = enable;
}