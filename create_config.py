import ConfigParser

config = ConfigParser.RawConfigParser()

# For the BRIGHT run
config.add_section('BRIGHT')
config.set('BRIGHT', '-c', 'config.sex')
config.set('BRIGHT', '-DETECT_MINAREA', '140')
config.set('BRIGHT', '-DETECT_THRESH', '2.2')
config.set('BRIGHT', '-DEBLEND_NTHRESH', '64')
config.set('BRIGHT', '-BACK_SIZE', '400')
config.set('BRIGHT', '-BACK_FILTERSIZE', '5')
config.set('BRIGHT', '-DEBLEND_MINCONT', '0.04')

# For the FAINT run
config.add_section('FAINT')
config.set('FAINT', '-c', 'config.sex')
config.set('FAINT', '-DETECT_MINAREA', '18')
config.set('FAINT', '-DETECT_THRESH', '1.0')
config.set('FAINT', '-DEBLEND_NTHRESH', '64')
config.set('FAINT', '-BACK_SIZE', '100')
config.set('FAINT', '-BACK_FILTERSIZE', '3')
config.set('FAINT', '-DEBLEND_MINCONT', '0.065')

# For the SMOOTHED/FAINT run
config.add_section('SMOOTH')
config.set('SMOOTH', '-c', 'config.sex')
config.set('SMOOTH', '-DETECT_MINAREA', '18')
config.set('SMOOTH', '-DETECT_THRESH', '1.0')
config.set('SMOOTH', '-DEBLEND_NTHRESH', '64')
config.set('SMOOTH', '-BACK_SIZE', '100')
config.set('SMOOTH', '-BACK_FILTERSIZE', '3')
config.set('SMOOTH', '-DEBLEND_MINCONT', '0.065')
config.set('SMOOTH', '-FILTER_NAME', 'sexfiles/filters/gauss_5.0_9x9.conv')


with open('se_params.cfg', 'wb') as configfile:
    config.write(configfile)


