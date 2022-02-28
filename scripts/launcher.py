
# internal function
def fire(appPath, graphPath, app, graph, version, policy):
    "just avoiding copy pasting"

    outFile = '../raw_data/out_' + app + '_' + graph + '_' + policy + '_' + version + '.dat'
    logFile = '../raw_data/log_' + app + '_' + graph + '_' + policy + '_' + version + '.log'

    runCmd = appPath + '/' + app + ' -f ' + graphPath + '/' + graph + '.sg ' \
             + ' -n 1 -i 1 2> ' + logFile + ' | tee  ' + outFile

    
    pinHeader = '../pin-3.21/pin -ifeellucky -t ../simulators/' + policy + '/cache_pinsim.so -- ' 
    
    runCmd = pinHeader + runCmd

    return runCmd


def launchSim(app, graph, policy, version):
    "launches run for an app after selecting \
    app specific cli options"
    
    graphPath = '../input-graphs/'
    
    srcPath = '../applications/' + version + '/'
    
    return fire(srcPath, graphPath, app, graph, version, policy)
    
