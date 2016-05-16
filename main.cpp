#include "MapBuilder.h"
#include <iostream>

int main()
{
    MapBuilder m(10,true,CERESBA);
    while(true){

        if(!m.process())
            break;
    }
}
