import logging
import sys


def create_logger(name,log_file,level=logging.INFO):

    try:
        handler=logging.FileHandler(log_file,mode='a+',encoding='utf-8')
        formatter=logging.Formatter('%(message)s')
        handler.setFormatter(formatter)
        logger=logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
    except PermissionError:
        sys.exit(-20)
    return logger

def log_to_user(logger,information, print_on_screen=1):
    if print_on_screen:
        print(information)
    logger.info(information)


def log_get_error(information,log_path):
    print(information)
    logging.basicConfig(filename=log_path,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.error(information)
    sys.exit(-10)


def log_get_warning(information,log_path):
    print(information)
    logging.basicConfig(filename=log_path,
                        filemode='a',
                        format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s'
                        )
    logging.warning(information)
