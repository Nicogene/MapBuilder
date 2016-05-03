#include "MapBuilder.h"
#include <iostream>

int main()
{
    MapBuilder m(5,true,CERESBA);
    while(true){

        if(!m.process())
            break;
    }
}
