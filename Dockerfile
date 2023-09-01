# Docker file that installs docker container for Selenzyme
#
# build with: "sudo docker build -t -i selenzyme_new ."
# sbc/selenzyme is the base image
FROM sbc/selenzybase2023

#ENV VIRTUAL_ENV=/opt/venv
#RUN python3 -m venv $VIRTUAL_ENV
#ENV PATH="$VIRTUAL_ENV/bin:$PATH"


ENTRYPOINT ["python"]
CMD ["/selenzyme2/selenzyPro/flaskform.py", "-uploaddir", "/selenzyme2/selenzyPro/uploads", "-datadir", "/selenzyme2/selenzyPro/data", "-logdir", "/selenzyme2/selenzyPro/log", "-d"] 

EXPOSE 5001

