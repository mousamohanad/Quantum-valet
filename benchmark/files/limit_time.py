import signal

def long_function_call():
    i = 0
    while i<1 :
	print(i)
        j=1
def signal_handler(signum, frame):
    raise Exception("Timed out!")

signal.signal(signal.SIGALRM, signal_handler)
signal.alarm(10)   # Ten seconds

try:
    long_function_call()
except Exception, msg:
    print "Timed out2!"
