#!bin/sh

find /home/ruth/code/update_selenzyme/selenzyme_2023/selenzyme2/selenzyPro/uploads/ -type f -mmin +15 | xargs rm

#rm -rf /home/ruth/code/update_selenzyme/selenzyme_2023/selenzyme2/selenzyPro/uploads/*
#find /home/ruth/code/update_selenzyme/selenzyme_2023/selenzyme2/selenzyPro/uploads/ -type f -mmin +1 -exec rm {} +
#find /home/ruth/code/update_selenzyme/selenzyme_2023/selenzyme2/selenzyPro/uploads/ -type f -mtime +1 -exec rm {} +
