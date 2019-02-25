import time

def printstatus(msg):
    '''
    Prints current step in the workflow.
    '''
    now = time.asctime( time.localtime(time.time()) )
    print("[{}]\t{}".format(now, msg))

def reportProgress(current,total):
    return "Completed: {0}% ({1} out of {2})".format( str(round( (current / total) * 100, 2)), current, total)
