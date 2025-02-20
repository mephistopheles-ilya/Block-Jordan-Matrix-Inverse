#! /bin/bash

falg=0

for((n=1;n<=10;++n))
do
    for((m=1;m<=n;++m))
    do
        for((s=1;s<=4;++s))
        do
            ../Thread/a.out ${n} ${m} 1 10 ${s} > lin_res.txt
            for((p=1;p<=6;++p))
            do
                mpirun -np ${p} -x ASAN_OPTIONS=detect_leaks=0 ./a.out ${n} ${m} 10 ${s} > mpi_res.txt
                diff mpi_res.txt lin_res.txt
                if [ $? -ne 0 ]; then
                    echo "Test Faild"
                    flag=1
                    echo ${n} ${m} 10 ${s}
                    break
                fi
            done
            if [[ $flag -ne 0 ]]; then
                break;
            fi
        done
        if [[ $flag -ne 0 ]]; then
            break;
        fi
    done
    if [[ $flag -ne 0 ]]; then
        break;
    fi
done



