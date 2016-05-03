#include "MapBuilder.h"
#include <iostream>

int main()
{
    MapBuilder m(5,true);
    while(true){

        if(!m.process())
            break;
    }
}
