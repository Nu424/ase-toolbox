from datetime import datetime
import logging


class ConditionalLogger:
    """
    ログ出力を条件によって制御するラッパークラス。

    enabled=False の場合、すべてのログ出力を無効化する。
    try-delta-g.py から移植。
    """

    def __init__(self, base_logger, enabled: bool = True):
        self.base_logger = base_logger
        self.enabled = enabled

    def __getattr__(self, name):
        """
        logger.info(), logger.warning() などの呼び出しを透過的に処理する。
        enabled=False の場合は何もしない関数を返す。
        """
        if not self.enabled:
            # ログが無効な場合は何もしない関数を返す
            return lambda *args, **kwargs: None

        # ログが有効な場合は元のloggerのメソッドを返す
        return getattr(self.base_logger, name)


def setup_logger(prefix: str = "calc"):
    """デバッグ用ログの設定"""
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"{prefix}_{timestamp}.log"

    logging.basicConfig(
        level=logging.DEBUG,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[
            logging.FileHandler(log_filename, encoding="utf-8"),
            logging.StreamHandler(),  # コンソールにも出力
        ],
    )

    base_logger = logging.getLogger(__name__)
    base_logger.info(f"デバッグログファイル: {log_filename}")
    return ConditionalLogger(base_logger, enabled=True)
